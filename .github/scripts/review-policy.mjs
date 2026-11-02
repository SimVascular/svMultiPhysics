import { pathToFileURL } from 'node:url';

const REPOSITORY = 'SimVascular/svMultiPhysics';
const TEAM = 'refactoring-team';
const TEAM_ID = 20019155;
const LABEL = 'OOP Refactor';
const LABEL_ID = 9482575472;
const APP_ID = 5288861;
const CHECK = 'review-policy';

export class GitHub {
  constructor(token, fetcher = fetch) {
    if (!token) throw new Error('GH_TOKEN is required.');
    this.token = token;
    this.fetcher = fetcher;
  }

  async request(path, { method = 'GET', body } = {}) {
    if (!path.startsWith('/') || path.startsWith('//')) throw new Error('Expected a relative API path.');
    const url = new URL(path, 'https://api.github.com');
    if (url.origin !== 'https://api.github.com') throw new Error('Unexpected API path origin.');
    const response = await this.fetcher(url, {
      method, redirect: 'error', signal: AbortSignal.timeout(30_000),
      headers: {
        Accept: 'application/vnd.github+json', Authorization: `Bearer ${this.token}`,
        'X-GitHub-Api-Version': '2022-11-28', 'Content-Type': 'application/json',
      },
      body: body === undefined ? undefined : JSON.stringify(body),
    });
    if (!response.ok) {
      const error = new Error(`GitHub API ${response.status}: ${method} ${url.pathname}`);
      error.status = response.status;
      throw error;
    }
    return response.status === 204 ? null : response.json();
  }

  async list(path, field) {
    const items = [];
    for (let page = 1; ; page++) {
      const separator = path.includes('?') ? '&' : '?';
      const response = await this.request(`${path}${separator}per_page=100&page=${page}`);
      const batch = field ? response[field] : response;
      if (!Array.isArray(batch)) throw new Error(`Expected an API list: ${path}`);
      items.push(...batch);
      if (batch.length < 100) return items;
    }
  }
}

const repoPath = `/repos/${REPOSITORY}`;
const teamPath = `/orgs/SimVascular/teams/${TEAM}`;
const openPulls = `${repoPath}/pulls?state=open&base=main`;

function prFingerprint(pr) {
  return JSON.stringify([pr.number, pr.user.id, pr.state, pr.head.sha, pr.base.ref,
    pr.draft, pr.updated_at, pr.labels.map(label => [label.id, label.name]).sort()]);
}

function reviewFingerprint(reviews) {
  return JSON.stringify(opinions(reviews).map(review => [review.id, review.user.id, review.state]));
}

function changesPolicy(files) {
  // The files endpoint truncates at 3,000 files; an incomplete list cannot grant an exception.
  return files.length >= 3000 || files.some(file => [file.filename, file.previous_filename]
    .some(path => path?.startsWith('.github/workflows/')
      || path?.startsWith('.github/actions/')
      || path?.startsWith('.github/scripts/review-policy')));
}

async function publishCheck(api, sha, passed, summary, detailsUrl) {
  let previous;
  try {
    const existing = await api.list(`${repoPath}/commits/${sha}/check-runs?check_name=${CHECK}&filter=latest`, 'check_runs');
    previous = existing.filter(check => check.name === CHECK && check.app?.id === APP_ID)
      .sort((a, b) => b.id - a.id)[0];
  } catch (error) {
    // A lookup failure must not prevent publishing a new failing check.
    if (passed) throw error;
  }
  const body = {
    name: CHECK, status: 'completed', conclusion: passed ? 'success' : 'failure',
    output: { title: passed ? 'Review policy satisfied' : 'Review policy not satisfied', summary },
    ...(detailsUrl ? { details_url: detailsUrl } : {}),
  };
  // Scheduled refreshes have different details URLs even when the result is unchanged.
  if (previous?.status === 'completed' && previous.conclusion === body.conclusion
      && previous.output?.title === body.output.title && previous.output?.summary === summary) return;
  if (previous && previous.status !== 'completed') {
    await api.request(`${repoPath}/check-runs/${previous.id}`, { method: 'PATCH', body });
  } else {
    await api.request(`${repoPath}/check-runs`, { method: 'POST', body: { ...body, head_sha: sha } });
  }
}

export async function reconcile(api, { apply = false, detailsUrl } = {}) {
  const prs = await api.list(openPulls);
  const heads = [...new Set(prs.map(pr => pr.head.sha))];
  if (heads.some(sha => !/^[a-f0-9]{40}$/.test(sha))) throw new Error('Invalid PR head SHA.');
  const records = [];
  try {
    // Reconcile every open PR in one serialized run: dropped queued events lose no PRs.
    if (!prs.length) return [];
    const configuredTeam = await api.request(teamPath);
    if (configuredTeam.id !== TEAM_ID) throw new Error('Refactoring team identity does not match configuration.');
    const team = new Set((await api.list(`${teamPath}/members`)).map(member => member.id));
    const permissions = new Map();
    for (const original of prs) {
      const path = `${repoPath}/pulls/${original.number}`;
      let pr = await api.request(path);
      const reviews = await api.list(`${path}/reviews`);
      const labeled = pr.labels.some(label => label.name === LABEL);
      const timeline = labeled ? await api.list(`${repoPath}/issues/${pr.number}/timeline`) : [];
      const files = labeled ? await api.list(`${path}/files`) : [];
      const writers = new Set();
      for (const { user } of opinions(reviews)) {
        if (user.id === pr.user.id) continue;
        if (!permissions.has(user.id)) {
          let eligible = false;
          try {
            const permission = await api.request(`${repoPath}/collaborators/${encodeURIComponent(user.login)}/permission`);
            eligible = permission.user?.id === user.id && ['write', 'admin'].includes(permission.permission);
          } catch (error) {
            if (error.status !== 404) throw error;
          }
          permissions.set(user.id, eligible);
        }
        if (permissions.get(user.id)) writers.add(user.id);
      }
      const result = { number: pr.number, head_sha: original.head.sha,
        ...evaluate({ pr, reviews, timeline, policyChanged: changesPolicy(files) }, team, writers) };
      if (apply && result.removeLabel) {
        const currentTimeline = await api.list(`${repoPath}/issues/${pr.number}/timeline`);
        const current = await api.request(path);
        if (prFingerprint(current) !== prFingerprint(pr) || JSON.stringify(currentTimeline) !== JSON.stringify(timeline)) {
          throw new Error(`PR #${pr.number} changed before label removal; run the policy again.`);
        }
        await api.request(`${repoPath}/issues/${pr.number}/labels/${encodeURIComponent(LABEL)}`, { method: 'DELETE' });
        pr = await api.request(path);
        result.reason += ' Unauthorized OOP Refactor label removed.';
      }
      if (pr.state !== 'open' || pr.base.ref !== 'main' || pr.head.sha !== original.head.sha) {
        result.passed = false;
        result.reason = 'PR changed during evaluation; run the policy again.';
      }
      records.push({ result, pr, reviews });
    }
    if (apply) {
      for (const record of records.filter(record => record.result.passed)) {
        const path = `${repoPath}/pulls/${record.pr.number}`;
        const reviews = await api.list(`${path}/reviews`);
        const pr = await api.request(path);
        if (prFingerprint(pr) !== prFingerprint(record.pr)
            || reviewFingerprint(reviews) !== reviewFingerprint(record.reviews)) {
          record.result.passed = false;
          record.result.reason = 'PR or reviews changed during evaluation; run the policy again.';
        }
      }
      const current = await api.list(openPulls);
      for (const sha of heads) {
        const group = records.filter(record => record.result.head_sha === sha).map(record => record.result)
          .sort((a, b) => a.number - b.number);
        const numbers = pulls => pulls.filter(pr => pr.head.sha === sha).map(pr => pr.number).sort((a, b) => a - b).join(',');
        const unchanged = numbers(current) === numbers(prs);
        const summary = group.map(result => `PR #${result.number}: ${result.reason}`).join('\n\n')
          + (unchanged ? '' : '\n\nPRs sharing this commit changed; run the policy again.');
        await publishCheck(api, sha, unchanged && group.every(result => result.passed), summary, detailsUrl);
      }
    }
    return records.map(record => record.result);
  } catch (error) {
    // Replace earlier successes too, including failures after one group was published.
    if (apply) await Promise.allSettled(heads.map(sha =>
      publishCheck(api, sha, false, `Evaluation could not complete. ${error.message}`, detailsUrl)));
    throw error;
  }
}

function opinions(reviews) {
  const latest = new Map();
  const submitted = reviews.filter(review => review.submitted_at);
  for (const review of submitted.sort((a, b) =>
    a.submitted_at.localeCompare(b.submitted_at) || a.id - b.id)) {
    if (review.user?.type !== 'User') continue;
    if (['APPROVED', 'CHANGES_REQUESTED', 'DISMISSED'].includes(review.state)) {
      latest.set(review.user.id, review);
    }
  }
  return [...latest.values()];
}

export function evaluate({ pr, reviews, timeline, policyChanged }, team, writers) {
  const label = pr.labels.find(item => item.name === LABEL);
  const application = [...timeline].reverse().find(event =>
    event.label?.name === LABEL && ['labeled', 'unlabeled'].includes(event.event));
  const authorized = label?.id === LABEL_ID && application?.event === 'labeled'
    && application.actor?.type === 'User' && team.has(application.actor.id);
  const exception = authorized && !policyChanged;
  const required = exception ? 1 : 2;
  const eligible = opinions(reviews).filter(review =>
    review.user.id !== pr.user.id && writers.has(review.user.id));
  const approvals = eligible.filter(review => review.state === 'APPROVED'
    && (!exception || team.has(review.user.id))).map(review => review.user.login);
  const blocked = eligible.some(review => review.state === 'CHANGES_REQUESTED');
  return {
    passed: !blocked && approvals.length >= required,
    required,
    approvals,
    removeLabel: Boolean(label && !authorized),
    reason: blocked ? 'An eligible reviewer has requested changes.'
      : `${approvals.length}/${required} ${exception ? 'Refactoring Team' : 'eligible'} approvals.`
        + (policyChanged ? ' Policy changes require two approvals.' : ''),
  };
}

async function main() {
  const args = process.argv.slice(2);
  if (args.some(arg => arg !== '--apply')) throw new Error('Usage: node review-policy.mjs [--apply]');
  const apply = args.includes('--apply');
  if (apply && (process.env.GITHUB_REPOSITORY !== REPOSITORY || process.env.GITHUB_REF !== 'refs/heads/main')) {
    throw new Error('Apply mode must run from SimVascular/svMultiPhysics main.');
  }
  const detailsUrl = process.env.GITHUB_RUN_ID
    ? `https://github.com/${REPOSITORY}/actions/runs/${process.env.GITHUB_RUN_ID}` : undefined;
  const results = await reconcile(new GitHub(process.env.GH_TOKEN), { apply, detailsUrl });
  console.log(JSON.stringify({ mode: apply ? 'apply' : 'dry-run', results }, null, 2));
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  main().catch(error => { console.error(error.message); process.exitCode = 1; });
}
