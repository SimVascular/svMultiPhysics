import assert from 'node:assert/strict';
import test from 'node:test';
import { evaluate, reconcile } from './review-policy.mjs';

const repo = '/repos/SimVascular/svMultiPhysics';
const teamPath = '/orgs/SimVascular/teams/refactoring-team';
const definitions = {
  'OOP Refactor': { name: 'OOP Refactor', id: 9482575472 },
  'software engineering': { name: 'software engineering', id: 12945018775 },
};
const user = id => ({ id, login: 'user-' + id, type: 'User' });
const review = (id, state = 'APPROVED') => ({
  id, user: user(id), state, submitted_at: '2026-11-03T00:00:00Z',
});
function fixture(names = ['software engineering'], actor = 1) {
  const labels = names.map(name => ({ ...definitions[name] }));
  return {
    pr: { number: 123, user: user(1), state: 'open', base: { ref: 'main' },
      head: { sha: 'a'.repeat(40) }, draft: false, updated_at: '2026-11-03T00:00:00Z', labels },
    reviews: [review(2)],
    timeline: labels.map(label => ({ event: 'labeled', label: { ...label }, actor: user(actor) })),
    policyChanged: false, filesTruncated: false,
    files: [{ filename: 'Code/example.cpp' }],
  };
}
function evaluateFixture(f) {
  return evaluate(f, new Set([1, 2]), new Set([1, 2, 3, 4]));
}
function fakeApi(f, afterDelete) {
  const deleted = [], checks = [];
  return {
    deleted, checks,
    async list(path) {
      if (path === repo + '/pulls?state=open&base=main') return [structuredClone(f.pr)];
      if (path === teamPath + '/members') return [user(1), user(2)];
      if (path === repo + '/pulls/123/reviews') return structuredClone(f.reviews);
      if (path === repo + '/pulls/123/files') return structuredClone(f.files);
      if (path === repo + '/issues/123/timeline') return structuredClone(f.timeline);
      if (path.includes('/check-runs?')) return [];
      throw new Error('Unexpected list: ' + path);
    },
    async request(path, options = {}) {
      if (path === teamPath) return { id: 20019155 };
      if (path === repo + '/pulls/123') return structuredClone(f.pr);
      if (path.startsWith(repo + '/collaborators/')) {
        const id = Number(path.match(/user-(\d+)/)[1]);
        return { permission: 'write', user: user(id) };
      }
      if (path.startsWith(repo + '/issues/123/labels/') && options.method === 'DELETE') {
        const name = decodeURIComponent(path.split('/').at(-1));
        deleted.push(name);
        f.pr.labels = f.pr.labels.filter(label => label.name !== name);
        f.timeline.push({ event: 'unlabeled', label: { ...definitions[name] },
          actor: { id: 99, type: 'Bot' } });
        afterDelete?.(f, name);
        return null;
      }
      if (path === repo + '/check-runs' && options.method === 'POST') {
        checks.push(options.body);
        return { id: checks.length };
      }
      throw new Error('Unexpected request: ' + path);
    },
  };
}
const pathCases = [
  ['source', 'Code/example.cpp', 1, 1],
  ['workflow', '.github/workflows/fast-tests.yml', 2, 1],
  ['action', '.github/actions/test-ubuntu/action.yml', 2, 1],
  ['review-policy workflow', '.github/workflows/review-policy.yml', 2, 1],
  ['review-event workflow', '.github/workflows/review-policy-events.yml', 2, 1],
  ['review-policy script', '.github/scripts/review-policy.mjs', 2, 1],
];
for (const [description, filename, oopRequired, softwareRequired] of pathCases) {
  for (const [name, expected] of [
    ['OOP Refactor', oopRequired], ['software engineering', softwareRequired],
  ]) {
    test(name + ': ' + description, async () => {
      const f = fixture([name]);
      f.files = [{ filename }];
      const [result] = await reconcile(fakeApi(f));
      assert.equal(result.required, expected);
      assert.equal(result.passed, expected === 1);
      assert.deepEqual(result.removeLabels, []);
      if (expected === 1) assert.doesNotMatch(result.reason, /require two approvals/);
    });
  }
}
test('unlabeled PR needs two approvals', () => {
  const f = fixture([]);
  assert.equal(evaluateFixture(f).passed, false);
  f.reviews.push(review(3));
  assert.equal(evaluateFixture(f).passed, true);
});
test('both authorized labels allow policy changes with one approval', () => {
  const f = fixture(Object.keys(definitions));
  f.policyChanged = true;
  assert.equal(evaluateFixture(f).required, 1);
});
test('unauthorized software engineering cannot override OOP restriction', () => {
  const f = fixture(Object.keys(definitions));
  f.timeline[1].actor = user(3);
  f.policyChanged = true;
  const result = evaluateFixture(f);
  assert.equal(result.required, 2);
  assert.deepEqual(result.removeLabels, ['software engineering']);
});
test('wrong label identity does not grant the exception', () => {
  const f = fixture();
  f.pr.labels[0].id++;
  assert.equal(evaluateFixture(f).required, 2);
  assert.deepEqual(evaluateFixture(f).removeLabels, ['software engineering']);
});
test('missing timeline does not grant the exception', () => {
  const f = fixture();
  f.timeline = [];
  assert.equal(evaluateFixture(f).required, 2);
});
test('bot-applied label does not grant the exception', () => {
  const f = fixture();
  f.timeline[0].actor.type = 'Bot';
  assert.equal(evaluateFixture(f).required, 2);
});
test('latest label application must be authorized', () => {
  const f = fixture();
  f.timeline.push({ ...f.timeline[0], event: 'unlabeled' });
  f.timeline.push({ ...f.timeline[0], actor: user(3) });
  assert.equal(evaluateFixture(f).required, 2);
});
test('approval must come from the team', () => {
  const f = fixture();
  f.reviews = [review(3)];
  assert.equal(evaluateFixture(f).passed, false);
});
test('author approval does not count', () => {
  const f = fixture();
  f.reviews = [review(1)];
  assert.equal(evaluateFixture(f).passed, false);
});
test('reviewer must have write permission', () => {
  const f = fixture();
  assert.equal(evaluate(f, new Set([1, 2]), new Set()).passed, false);
});
test('changes requested still block an exception', () => {
  const f = fixture();
  f.reviews.push(review(3, 'CHANGES_REQUESTED'));
  assert.equal(evaluateFixture(f).passed, false);
});
test('renamed policy file still restricts OOP Refactor', async () => {
  const f = fixture(['OOP Refactor']);
  f.files = [{ filename: 'old-policy.mjs', previous_filename: '.github/scripts/review-policy.mjs' }];
  const [result] = await reconcile(fakeApi(f));
  assert.equal(result.required, 2);
});
test('truncated file list does not grant software engineering exception', async () => {
  const f = fixture();
  f.files = Array.from({ length: 3000 }, () => ({ filename: 'Code/example.cpp' }));
  const [result] = await reconcile(fakeApi(f));
  assert.equal(result.required, 2);
  assert.match(result.reason, /incomplete/);
});
test('dry run identifies both unauthorized labels without removing them', async () => {
  const f = fixture(Object.keys(definitions), 3);
  const api = fakeApi(f);
  const [result] = await reconcile(api);
  assert.deepEqual(result.removeLabels, Object.keys(definitions));
  assert.deepEqual(api.deleted, []);
  assert.deepEqual(api.checks, []);
});
test('apply removes both unauthorized labels and publishes the reevaluated result', async () => {
  const f = fixture(Object.keys(definitions), 3);
  f.reviews.push(review(3));
  const api = fakeApi(f);
  const [result] = await reconcile(api, { apply: true });
  assert.deepEqual(api.deleted, Object.keys(definitions));
  assert.deepEqual(result.removeLabels, []);
  assert.equal(result.required, 2);
  assert.equal(result.passed, true);
  assert.equal(api.checks[0].conclusion, 'success');
});
test('removing unauthorized OOP label preserves authorized software engineering exception', async () => {
  const f = fixture(Object.keys(definitions));
  f.timeline[0].actor = user(3);
  f.files = [{ filename: '.github/workflows/review-policy.yml' }];
  const api = fakeApi(f);
  const [result] = await reconcile(api, { apply: true });
  assert.deepEqual(api.deleted, ['OOP Refactor']);
  assert.equal(result.required, 1);
  assert.equal(result.passed, true);
  assert.equal(api.checks[0].conclusion, 'success');
});
test('label authorized between removals is reevaluated and preserved', async () => {
  const f = fixture(Object.keys(definitions), 3);
  f.files = [{ filename: '.github/scripts/review-policy.mjs' }];
  const api = fakeApi(f, (current, name) => {
    if (name === 'OOP Refactor') current.timeline.push({
      event: 'labeled', label: { ...definitions['software engineering'] }, actor: user(1),
    });
  });
  const [result] = await reconcile(api, { apply: true });
  assert.deepEqual(api.deleted, ['OOP Refactor']);
  assert.equal(result.required, 1);
  assert.equal(result.passed, true);
});
