# Cardiac Fiber Generation Documentation

This document describes the mathematical framework for generating myocardial fiber orientations in cardiac geometries using Laplace-Dirichlet rule-based methods.

# General Overview of the Fiber Generation Methods
To explain the main steps and functions of the code, in the following we explain the main concepts and methods to generate fibers. In the next section these concepts are used to describe the Bayer and Doste method. 

## 1. Laplace Problem

The foundation of the rule-based fiber generation is the solution of Laplace-Dirichlet boundary value problems. These provide scalar fields $\phi$ (surface1 → surface2) that satisfy:

$$
\begin{align}
\nabla^2 \phi = 0 &\quad \text{in } \Omega \\
\phi_t = 0 & \quad\text{on } S_{\text{surface1}} \\
\phi_t = 1 & \quad\text{on } S_{\text{surface2}}
\end{align}
$$

### 1.1 Transmural Direction

The transmural field $\phi_t$ (endo → epi) characterizes the wall thickness direction, varying from endocardium to epicardium:

$$
\begin{cases}
\nabla^2 \phi_t = 0 & \text{in } \Omega \\
\phi_t = 0 & \text{on } S_{\text{endo}} \\
\phi_t = 1 & \text{on } S_{\text{epi}}
\end{cases}
$$

This field is normalized to $[0, 1]$ range where $\phi_t = 0$ at the endocardium and $\phi_t = 1$ at the epicardium.

### 1.2 Longitudinal Direction

The longitudinal field $\phi_\ell$ (apex → base) characterizes the apex-to-base direction:

$$
\begin{cases}
\nabla^2 \phi_\ell = 0 & \text{in } \Omega \\
\phi_\ell = 0 & \text{on } S_{\text{apex}} \\
\phi_\ell = 1 & \text{on } S_{\text{base}}
\end{cases}
$$

This field is also normalized to $[0, 1]$ where $\phi_\ell = 0$ at the apex and $\phi_\ell = 1$ at the base.

**Implementation Note**: The Laplace equations are solved using the SVMultiphysics solver configured to solve a steady-state heat equation (which is equivalent to the Laplace equation). The solver is configured with:
- `Conductivity = 1.0`, `Source_term = 0.0`, `Density = 0.0`
- `Spectral_radius_of_infinite_time_step = 0.0`
- Single time step to obtain the steady-state solution directly

## 2. Definition of Basis

A local orthonormal basis $\{\mathbf{e}_c, \mathbf{e}_\ell, \mathbf{e}_t\}$ is constructed at each point in the myocardium, where:
- $\mathbf{e}_c$: circumferential direction
- $\mathbf{e}_\ell$: longitudinal direction  
- $\mathbf{e}_t$: transmural direction

### 2.1 Obtain Gradients from Laplace Solutions

The gradients of the Laplace fields provide natural directional vectors. The gradients are computed at mesh nodes and then averaged to cell centers for smoother results:

$$
\mathbf{g}_t = \nabla \phi_t, \quad \mathbf{g}_\ell = \nabla \phi_\ell
$$

These gradients are normalized to unit vectors:

$$
\hat{\mathbf{g}}_t = \frac{\mathbf{g}_t}{\|\mathbf{g}_t\|}, \quad \hat{\mathbf{g}}_\ell = \frac{\mathbf{g}_\ell}{\|\mathbf{g}_\ell\|}
$$

### 2.2 Calculate Circumferential Direction

The local orthonormal basis is constructed using the `axis` function. Given the longitudinal direction $\hat{\mathbf{g}}_\ell$ and the transmural direction $\hat{\mathbf{g}}_t$:

1. **Longitudinal basis vector**: 
    $$\mathbf{e}_\ell = \frac{\hat{\mathbf{g}}_\ell}{\|\hat{\mathbf{g}}_\ell\|}$$

2. **Transmural basis vector** (orthogonalized to $\mathbf{e}_\ell$):
    $$\mathbf{e}_t' = \hat{\mathbf{g}}_t - (\hat{\mathbf{g}}_t \cdot \mathbf{e}_\ell)\mathbf{e}_\ell$$
    $$\mathbf{e}_t = \frac{\mathbf{e}_t'}{\|\mathbf{e}_t'\|}$$

3. **Circumferential basis vector** (orthogonal to both):
    $$\mathbf{e}_c = \mathbf{e}_\ell \times \mathbf{e}_t$$

This ensures a right-handed orthonormal coordinate system at each element. We implement this in the **axis** function, $\mathbf Q=(\mathbf{e}_c, \mathbf{e}_\ell, \mathbf{e}_t)=\text{axis}(\hat{\mathbf{g}}_\ell, \hat{\mathbf{g}}_t)$.

## 3. Definition of Angles over the Geometry

Two angles define the fiber orientation relative to the local basis $\mathbf Q=(\mathbf e_c, \mathbf e_\ell, \mathbf e_t)$:

- **$\alpha$ (helix angle)**: Rotation angle using as axis of rotation $\mathbf{e}_t$. Returns a rotated basis $\mathbf Q^{(\alpha)}=(\mathbf e_c', \mathbf e_\ell', \mathbf e_t)$
- **$\beta$ (transverse angle)**: Rotation angle using as axis of rotation $\mathbf{e}_\ell'$. Returns the fiber ($\mathbf f$), sheetnormal ($\mathbf n$), and sheet ($\mathbf s$) vectors $\mathbf Q^{(\alpha, \beta)} = (\mathbf f, \mathbf n, \mathbf s)$.

The angles vary linearly across the wall thickness based on the transmural coordinate $\phi_t$. For the single ventricle this is:

$$
\alpha(\phi_t) = \alpha_{\text{endo}}(1 - \phi_t) + \alpha_{\text{epi}}\phi_t
$$

$$
\beta(\phi_t) = \beta_{\text{endo}}(1 - \phi_t) + \beta_{\text{epi}}\phi_t
$$

where:
- $\alpha_{\text{endo}}, \alpha_{\text{epi}}$: helix angles at endocardium and epicardium
- $\beta_{\text{endo}}, \beta_{\text{epi}}$: transverse angles at endocardium and epicardium


## 4. Rotation of the Basis

The fiber direction is obtained by applying two successive rotations to the local basis.

### 4.1 Rotation Using Matrices (Bayer)

The two-step rotation can be expressed with standard rotation matrices applied to the local basis $\mathbf{Q} = [\mathbf{e}_c, \mathbf{e}_\ell, \mathbf{e}_t]$:

- Helix rotation by $\alpha$ about the transmural axis $\mathbf{e}_t$:

$$
\mathbf{R}_\alpha = \begin{bmatrix}
\cos\alpha & -\sin\alpha & 0 \\
\sin\alpha & \cos\alpha & 0 \\
0 & 0 & 1
\end{bmatrix}, \quad
\mathbf{Q}^{(\alpha)} = \mathbf{Q}\,\mathbf{R}_\alpha.
$$

- Transverse rotation by $\beta$ about the rotated longitudinal axis $\mathbf{e}_\ell^{(\alpha)}$ (the second column of $\mathbf{Q}^{(\alpha)}$):

$$
\mathbf{R}_\beta = \begin{bmatrix}
\cos\beta & 0 & \sin\beta \\
0 & 1 & 0 \\
-\sin\beta & 0 & \cos\beta
\end{bmatrix}, \quad
\mathbf{Q}^{(\alpha,\beta)} = \mathbf{Q}^{(\alpha)}\,\mathbf{R}_\beta.
$$

We implement this in the **orient\_matrix** function, $\mathbf Q^{(\alpha, \beta)}=(\mathbf{f}, \mathbf{n}, \mathbf{s})=\text{orient\_matrix}(\mathbf Q, \alpha, \beta)$.

Note: In the original Bayer paper the second matrix was written as,
$$
\mathbf{R}_\beta = \begin{bmatrix}
1 & 0 & 0 \\
0 & \cos\beta &  \sin\beta \\
0 & -\sin\beta &  \cos\beta 
\end{bmatrix}.
$$
Given how the orthogonal basis are ordered $(\mathbf{e}_c', \mathbf{e}_\ell', \mathbf{e}_t)$, this is equivalent to rotate over the first vector $\mathbf{e}_c'$ which is not what we want and keeps the fiber orientation unchanged and unaffected by $\beta$ angles.


### 4.2 Rotation Using Rodrigues' Formula (Doste)

These rotations can also be achieved using the Rodrigues rotation formula. For a unit axis $\mathbf{n}$ and angle $\theta$:

$$
\mathbf{R}(\theta, \mathbf{n}) = \mathbf{I}\cos\theta + [\mathbf{n}]_\times \sin\theta + \mathbf{n}\,\mathbf{n}^T (1 - \cos\theta),
$$

where the skew-symmetric matrix $[\mathbf{n}]_\times$ is

$$
[\mathbf{n}]_\times = \begin{bmatrix}
0 & -n_z & n_y \\
n_z & 0 & -n_x \\
-n_y & n_x & 0
\end{bmatrix}, \quad \mathbf{n} = (n_x, n_y, n_z)^T.
$$

Applying the two-step rotation to the local basis:

- Helix rotation by $\alpha$ about $\mathbf{e}_t$:

$$
\mathbf{Q}^{(\alpha)} = \mathbf{Q}\,\mathbf{R}(\alpha, \mathbf{e}_t).
$$

- Transverse rotation by $\beta$ about the rotated longitudinal axis $\mathbf{e}_\ell^{(\alpha)}$:

$$
\mathbf{Q}^{(\alpha,\beta)} = \mathbf{Q}^{(\alpha)}\,\mathbf{R}(\beta, \mathbf{e}_\ell^{(\alpha)}), \quad \mathbf{e}_\ell^{(\alpha)} = \big(\mathbf{Q}^{(\alpha)}\big)[:, 1].
$$

We implement this in the **orient\_rodrigues** function, $\mathbf Q^{(\alpha, \beta)}=(\mathbf{f}, \mathbf{n}, \mathbf{s})=\text{orient\_rodrigues}(\mathbf Q, \alpha, \beta)$.

## 5. Basis Interpolation

When working with biventricular geometries, different orthogonal basis are computed for the left ventricle (LV) and right ventricle (RV). For this, the orthogonal basis are represented as quaternions, which are then interpolated using spherical linear interpolation (SLERP).

### Spherical Linear Interpolation (SLERP)

Simple linear interpolation of rotation matrices can produce non-orthogonal results. Instead, **interpolate_basis** (bilinear spherical interpolation) is used, which operates via quaternion representation:

Given two rotation matrices $\mathbf{Q}_1$ and $\mathbf{Q}_2$, and interpolation parameter $t \in [0, 1]$:

1. **Convert to quaternions**: 
   $$\mathbf{q}_1 = \text{rotm2quat}(\mathbf{Q}_1), \quad \mathbf{q}_2 = \text{rotm2quat}(\mathbf{Q}_2)$$

2. **Ensure shortest path** (quaternion double cover):
   $$\text{if } \mathbf{q}_1 \cdot \mathbf{q}_2 < 0: \quad \mathbf{q}_2 \leftarrow -\mathbf{q}_2$$

3. **SLERP formula**:
   $$\theta_0 = \arccos(\mathbf{q}_1 \cdot \mathbf{q}_2)$$
   $$\mathbf{q}(t) = \frac{\sin((1-t)\theta_0)}{\sin\theta_0}\mathbf{q}_1 + \frac{\sin(t\theta_0)}{\sin\theta_0}\mathbf{q}_2$$

4. **Convert back to rotation matrix**:
   $$\mathbf{Q}(t) = \text{quat2rotm}(\mathbf{q}(t))$$

For nearly parallel quaternions ($\sin\theta_0 < 10^{-6}$), linear interpolation is used instead to avoid numerical issues.

We implement this in the **interpolate_basis** function, $\mathbf{Q}(t) = \text{interpolate\_basis}(\mathbf{Q}_1, \mathbf{Q}_2, t)$.

# Bayer Method

The Bayer et al. (2012) method is designed for truncated biventricular geometries without outflow tracts.

## Required Boundaries

The fiber generation algorithm requires specific boundary surfaces to be defined on the cardiac mesh:

- **Epicardium**: The outer surface of the biventricle
- **LV Endocardium**: The inner surface of the left ventricle
- **RV Endocardium**: The inner surface of the right ventricle
- **Base**: The basal (top) boundary of the geometry
- **Apex**: The epicardial apex region

The apex surface is automatically generated from the epicardium by identifying the point furthest from the base. Specifically, the apex point $\mathbf{p}_{\text{apex}}$ is found as:

$$
\mathbf{p}_{\text{apex}} = \arg\min_{\mathbf{p} \in S_{\text{epi}} \setminus S_{\text{base}}} \|\mathbf{p} - \mathbf{c}_{\text{base}}\|
$$

where $S_{\text{epi}}$ is the epicardial surface, $S_{\text{base}}$ is the base surface, and $\mathbf{c}_{\text{base}}$ is the centroid of the base. The apex surface consists of all elements in the epicardium that contain this apex point.

### Required Laplace Fields

Four Laplace problems are solved:

1. **Epi transmural**: $\phi_{\text{epi}}$ (LV endo and RV endo → epi)
2. **LV transmural**: $\phi_{\text{LV}}$ (RV endo and epi → LV endo)
2. **RV transmural**: $\phi_{\text{RV}}$ (LV endo and epi → RV endo)
4. **Apex-to-base**: $\phi_{\text{AB}}$ (apex → base)

## Input Angles

The Bayer method requires four input angle parameters (typically specified in degrees):

- **$\alpha_{\text{endo}}$**: Endocardial helix angle, typically $60°$
- **$\alpha_{\text{epi}}$**: Epicardial helix angle, typically $-60°$
- **$\beta_{\text{endo}}$**: Endocardial transverse angle, typically $-20°$
- **$\beta_{\text{epi}}$**: Epicardial transverse angle, typically $20°$

These angles define the fiber architecture that varies smoothly from endocardium to epicardium across both ventricles.

### Angle Definition

The angles vary based on the interventricular coordinate $d$ and transmural coordinate $\phi_{\text{epi}}$:

$$
d = \frac{\phi_{\text{RV}}}{\phi_{\text{LV}} + \phi_{\text{RV}}}
$$

**Septum angles** (interpolated between LV and RV):
$$
\alpha_s = \alpha_{\text{endo}}(1 - d) - \alpha_{\text{endo}} d
$$
$$
\beta_s = \beta_{\text{endo}}(1 - d) - \beta_{\text{endo}} d
$$

**Wall angles** (transmural variation):
$$
\alpha_w = \alpha_{\text{endo}}(1 - \phi_{\text{epi}}) + \alpha_{\text{epi}}\phi_{\text{epi}}
$$
$$
\beta_w = \beta_{\text{endo}}(1 - \phi_{\text{epi}}) + \beta_{\text{epi}}\phi_{\text{epi}}
$$

### Algorithm Steps (as in Bayer's paper)

1. **Construct LV and RV basis**:
   - $\mathbf{Q}_{\text{LV}}^0 = \text{axis}(\nabla\phi_{\text{AB}}, -\nabla\phi_{\text{LV}})$
   - $\mathbf{Q}_{\text{RV}}^0 = \text{axis}(\nabla\phi_{\text{AB}}, \nabla\phi_{\text{RV}})$

2. **Rotate by septum angles**:
   - $\mathbf{Q}_{\text{LV}} = \text{orient\_matrix}(\mathbf{Q}_{\text{LV}}^0, \alpha_s, \beta_s)$
   - $\mathbf{Q}_{\text{RV}} = \text{orient\_matrix}(\mathbf{Q}_{\text{RV}}^0, \alpha_s, \beta_s)$

3. **Interpolate endocardial basis**:
   - $\mathbf{Q}_{\text{endo}} = \text{bislerp}(\mathbf{Q}_{\text{LV}}, \mathbf{Q}_{\text{RV}}, d)$

4. **Construct epicardial basis**:
   - $\mathbf{Q}_{\text{epi}}^0 = \text{axis}(\nabla\phi_{\text{AB}}, \nabla\phi_{\text{epi}})$

5. **Rotate by wall angles**:
   - $\mathbf{Q}_{\text{epi}} = \text{orient\_matrix}(\mathbf{Q}_{\text{epi}}^0, \alpha_w, \beta_w)$

6. **Interpolate final basis**:
   - $\mathbf{Q} = \text{bislerp}(\mathbf{Q}_{\text{endo}}, \mathbf{Q}_{\text{epi}}, \phi_{\text{epi}})$
   - Extract: $\mathbf{f} = \mathbf{Q}[:, 0]$, $\mathbf{n} = \mathbf{Q}[:, 1]$, $\mathbf{s} = \mathbf{Q}[:, 2]$


Note. The bislerp function performs the same SLERP operation as interpolate_basis, but includes an additional correction to account for the fact that fiber directions are equivalent up to sign; that is, for the physical applications considered, using $\mathbf{f}$ or $-\mathbf{f}$ is equivalent. This correction flips each vector of the basis $\mathbf{Q}_1$ and selects the flipped configuration that is closest to the target basis $\mathbf{Q}_2$. However, this approach can introduce issues (see the next section), particularly when $\mathbf{Q}_1$ and $\mathbf{Q}_2$ are nearly orthogonal. In such cases, small perturbations in the basis can lead to drastically different interpolation results.

### Modified algorithm steps
When running the original implementation, we observed the resulting fibers showed discontinuities that arise due to the **bislerp** function. To solve these issues, we modify the algorithm as follows:


   - In step 2, we do $\mathbf{Q}_{\text{LV}} = \text{orient\_matrix}(\mathbf{Q}_{\text{LV}}^0, \alpha_s, \text{abs}(\beta_s))$ and $\mathbf{Q}_{\text{RV}} = \text{orient\_matrix}(\mathbf{Q}_{\text{LV}}^0, \alpha_s, \text{abs}(\beta_s))$.

      Note that $\mathbf{Q}_{\mathrm{LV}}^{0}$ and $\mathbf{Q}_{\mathrm{RV}}^{0}$ share equivalent circumferential, longitudinal, and transmural directions \textbf{within the septum}. By definition, $\beta_s > 0$ on the LV side (assuming $\beta_{\mathrm{endo}} > 0$), causing the fiber vector to rotate outward from the septum. On the RV side, $\beta_s < 0$ which also causes the fiber vector to rotate away from the septum. However, this is a negative angle at the RV endocardium, which is not what we want. Taking the absolute value $|\beta_s|$ yields the correct fiber angles while preserving the transmural variation of $\beta_{\mathrm{endo}}$ (positive at both side of the septum and 0 at the center of the septum).


      ![Illustration of beta angle effect at endocardium](example/biv_truncated/betaendo.png)
   
   - After step 3, for elements where $d > 0.5$, flip the first and third basis vectors (fiber and sheet) of $\mathbf{Q}_{\text{endo}}$. 
   
      Note that $\mathbf{Q}_{\text{LV}}^{0}$ and $\mathbf{Q}_{\text{RV}}^{0}$ are constructed with opposite signs for the transmural direction. As a result, the LV basis rotates counterclockwise, whereas the RV basis rotates clockwise. At the septum, the circumferential vectors of both bases point in the same direction, which allows for a straightforward SLERP interpolation along the shortest path to obtain $\mathbf{Q}_{\text{endo}}$. However, on the RV side this construction causes the $\mathbf{Q}_{\text{endo}}$ basis to point exactly opposite to $\mathbf{Q}_{\text{epi}}$, leading to issues with the SLERP interpolation. Flipping the vectors on the RV side resolves this problem and ensures that the second interpolation remains smooth.

            
      ![Illustration of beta angle effect at endocardium](example/biv_truncated/flipping.png)


# Doste Method

The Doste et al. (2019) method extends the fiber generation to biventricular geometries with **outflow tracts**, requiring valve surfaces to be explicitly defined.

## Required Boundaries

The fiber generation algorithm requires the following boundary surfaces:

- **Epicardium**: The outer surface of the biventricle
- **LV Endocardium**: The inner surface of the left ventricle
- **RV Endocardium**: The inner surface of the right ventricle
- **Mitral valve**: LV inflow boundary
- **Aortic valve**: LV outflow boundary
- **Tricuspid valve**: RV inflow boundary
- **Pulmonary valve**: RV outflow boundary
- **Apex**: The epicardial apex region
- **Top**: Surface that includes all valves. Only needed if the Apex needs to be generated.

The apex surface is automatically generated from the epicardium by identifying the point furthest from the base/top boundary. The apex point $\mathbf{p}_{\text{apex}}$ is found as:

$$
\mathbf{p}_{\text{apex}} = \arg\min_{\mathbf{p} \in S_{\text{epi}} \setminus S_{\text{top}}} \|\mathbf{p} - \mathbf{c}_{\text{top}}\|
$$

where $S_{\text{epi}}$ is the epicardial surface, $S_{\text{top}}$ is the top surface, and $\mathbf{c}_{\text{top}}$ is the centroid of the top boundary.

### Required Laplace Fields

Eleven Laplace problems are solved:

**Interventricular**:
1. $\phi_{\text{BiV}}$: RV endocardium → LV endocardium (interventricular field)
2. $\phi_{\text{Trans}}$: Combined transmural field for septal angle calculation

**Left Ventricle**:
1. $\phi_{\text{LV,trans}}$: Epicardium → LV endocardium (LV transmural for basis construction)
2. $\phi_{\text{LV,av}}$: Aortic valve → apex (LV longitudinal from AV)
3. $\phi_{\text{LV,mv}}$: Mitral valve → apex (LV longitudinal from MV)
4. $\phi_{\text{LV,weight}}$: Aortic valve → mitral valve (LV valve weight)

**Right Ventricle**:
1. $\phi_{\text{RV,trans}}$: Epicardium → RV endocardium (RV transmural for basis construction)
2. $\phi_{\text{RV,pv}}$: Pulmonary valve → apex (RV longitudinal from PV)
3. $\phi_{\text{RV,tv}}$: Tricuspid valve → apex (RV longitudinal from TV)
4. $\phi_{\text{RV,weight}}$: Pulmonary valve → tricuspid valve (RV valve weight)

**Global Transmural**:
1. $\phi_{\text{epi,trans}}$: LV and RV endocardium → epicardium (global transmural for angle interpolation)

**Note**: All Laplace fields are automatically normalized to the $[0, 1]$ range **except** for $\phi_{\text{BiV}}$ and $\phi_{\text{Trans}}$, which are kept in their original range for proper septal field calculation and basis construction. The $\phi_{\text{BiV}}$ field typically has negative values on the RV side and positive values on the LV side.

## Input Angles

The Doste method requires twelve input angle parameters (typically specified in degrees) to handle both ventricles and outflow tracts:

- **$\alpha_{\text{endo,LV}}$**: LV endocardial helix angle, typically $60°$
- **$\alpha_{\text{epi,LV}}$**: LV epicardial helix angle, typically $-60°$
- **$\beta_{\text{endo,LV}}$**: LV endocardial transverse angle, typically $-20°$
- **$\beta_{\text{epi,LV}}$**: LV epicardial transverse angle, typically $20°$
- **$\alpha_{\text{endo,RV}}$**: RV endocardial helix angle, typically $90°$
- **$\alpha_{\text{epi,RV}}$**: RV epicardial helix angle, typically $-25°$
- **$\beta_{\text{endo,RV}}$**: RV endocardial transverse angle, typically $-20°$
- **$\beta_{\text{epi,RV}}$**: RV epicardial transverse angle, typically $20°$
- **$\alpha_{\text{OT,endo,LV}}$**: LV outflow tract endocardial helix angle, typically $90°$
- **$\alpha_{\text{OT,epi,LV}}$**: LV outflow tract epicardial helix angle, typically $0°$
- **$\alpha_{\text{OT,endo,RV}}$**: RV outflow tract endocardial helix angle, typically $90°$
- **$\alpha_{\text{OT,epi,RV}}$**: RV outflow tract epicardial helix angle, typically $0°$

The outflow tract angles are blended with ventricular angles using the valve weight functions to create smooth transitions between the ventricular and outflow tract regions.

### Valve Weight Functions

To localize the influence of each valve, weight functions are computed and redistributed:

$$
w_{\text{LV}} = \text{redistribute}(\phi_{\text{LV,mv}}, q_{\text{up}} = 0.7, q_{\text{low}} = 0.01)
$$
$$
w_{\text{RV}} = \text{redistribute}(\phi_{\text{RV,tv}}, q_{\text{up}} = 0.1, q_{\text{low}} = 0.001)
$$

The redistribution clips values at quantiles $q_{\text{low}}$ and $q_{\text{up}}$, then renormalizes to $[0, 1]$.

### Angle Definition

The angles are region-specific with outflow tract blending:

**LV wall angles**:
$$
\alpha_{\text{LV,endo}} = \alpha_{\text{endo,LV}} w_{\text{LV}} + \alpha_{\text{OT,endo,LV}}(1 - w_{\text{LV}})
$$
$$
\alpha_{\text{LV,epi}} = \alpha_{\text{epi,LV}} w_{\text{LV}} + \alpha_{\text{OT,epi,LV}}(1 - w_{\text{LV}})
$$
$$
\alpha_{\text{wall,LV}} = \alpha_{\text{LV,endo}}(1 - \phi_{\text{epi,trans}}) + \alpha_{\text{LV,epi}}\phi_{\text{epi,trans}}
$$
$$
\beta_{\text{wall,LV}} = \left[\beta_{\text{endo,LV}}(1 - \phi_{\text{epi,trans}}) + \beta_{\text{epi,LV}}\phi_{\text{epi,trans}}\right] w_{\text{LV}}
$$

**RV wall angles**:
$$
\alpha_{\text{RV,endo}} = \alpha_{\text{endo,RV}} w_{\text{RV}} + \alpha_{\text{OT,endo,RV}}(1 - w_{\text{RV}})
$$
$$
\alpha_{\text{RV,epi}} = \alpha_{\text{epi,RV}} w_{\text{RV}} + \alpha_{\text{OT,epi,RV}}(1 - w_{\text{RV}})
$$
$$
\alpha_{\text{wall,RV}} = \alpha_{\text{RV,endo}}(1 - \phi_{\text{epi,trans}}) + \alpha_{\text{RV,epi}}\phi_{\text{epi,trans}}
$$
$$
\beta_{\text{wall,RV}} = \left[\beta_{\text{endo,RV}}(1 - \phi_{\text{epi,trans}}) + \beta_{\text{epi,RV}}\phi_{\text{epi,trans}}\right] w_{\text{RV}}
$$

**Septum angles**: Computed by blending LV and RV contributions weighted by septal position.

First, compute a septal field $s$ that is 1 at both endocardia but assigns 2/3 of the septum to the LV:
$$
s = \begin{cases}
\phi_{\text{Trans}} / 2 & \text{if } \phi_{\text{Trans}} < 0 \\
\phi_{\text{Trans}} & \text{otherwise}
\end{cases}
$$
$$
s = |s|
$$

Then compute septum angles based on which ventricle the element belongs to:
$$
\alpha_{\text{septum}} = \begin{cases}
\alpha_{\text{LV,endo}} \cdot s & \text{if } \phi_{\text{BiV}} < 0 \\
\alpha_{\text{RV,endo}} \cdot s & \text{if } \phi_{\text{BiV}} > 0
\end{cases}
$$
$$
\beta_{\text{septum}} = \begin{cases}
\beta_{\text{endo,LV}} \cdot s \cdot w_{\text{LV}} & \text{if } \phi_{\text{BiV}} < 0 \\
\beta_{\text{endo,RV}} \cdot s \cdot w_{\text{RV}} & \text{if } \phi_{\text{BiV}} > 0
\end{cases}
$$

### Algorithm Steps

1. **Construct LV and RV longitudinal directions**:
   - Blend valve gradients: $\mathbf{g}_{\ell,\text{LV}} = w_{\text{LV}} \nabla\phi_{\text{LV,mv}} + (1 - w_{\text{LV}})\nabla\phi_{\text{LV,av}}$
   - Blend valve gradients: $\mathbf{g}_{\ell,\text{RV}} = w_{\text{RV}} \nabla\phi_{\text{RV,tv}} + (1 - w_{\text{RV}})\nabla\phi_{\text{RV,pv}}$

2. **Construct LV and RV basis** using chamber-specific transmural fields:
   - $\mathbf{Q}_{\text{LV}} = \text{axis}(-\mathbf{g}_{\ell,\text{LV}}, -\nabla\phi_{\text{LV,trans}})$
   - $\mathbf{Q}_{\text{RV}} = \text{axis}(-\mathbf{g}_{\ell,\text{RV}}, -\nabla\phi_{\text{RV,trans}})$
   
   Note: The minus signs ensure the basis vectors follow the FibGen convention (apex → base, endo → epi).

3. **Rotate by septum angles to create septal bases**:
   - $\mathbf{Q}_{\text{LV,septum}} = \text{orient\_rodrigues}(\mathbf{Q}_{\text{LV}}, \alpha_{\text{septum}}, \beta_{\text{septum}})$
   - $\mathbf{Q}_{\text{RV,septum}} = \text{orient\_rodrigues}(\mathbf{Q}_{\text{RV}}, \alpha_{\text{septum}}, \beta_{\text{septum}})$

4. **Rotate by wall angles to create wall bases**:
   - $\mathbf{Q}_{\text{LV,wall}} = \text{orient\_rodrigues}(\mathbf{Q}_{\text{LV}}, \alpha_{\text{wall,LV}}, \beta_{\text{wall,LV}})$
   - $\mathbf{Q}_{\text{RV,wall}} = \text{orient\_rodrigues}(\mathbf{Q}_{\text{RV}}, \alpha_{\text{wall,RV}}, \beta_{\text{wall,RV}})$

5. **Create discontinuous septal basis**:
   
   Normalize $\phi_{\text{BiV}}$ for interpolation:
   $$
   \phi_{\text{BiV}}^{\text{norm}} = \begin{cases}
   \phi_{\text{BiV}} / 2 & \text{if } \phi_{\text{BiV}} < 0 \\
   \phi_{\text{BiV}} & \text{otherwise}
   \end{cases}
   $$
   $$
   \phi_{\text{BiV}}^{\text{interp}} = \frac{\phi_{\text{BiV}}^{\text{norm}} + 1}{2}
   $$
   
   Then assign septal basis discontinuously:
   $$\mathbf{Q}_{\text{sep}} = \begin{cases}
   \mathbf{Q}_{\text{LV,septum}} & \text{if } \phi_{\text{BiV}}^{\text{interp}} < 0.5 \\
   \mathbf{Q}_{\text{RV,septum}} & \text{if } \phi_{\text{BiV}}^{\text{interp}} \geq 0.5
   \end{cases}$$

6. **Interpolate epicardial basis**:
   - $\mathbf{Q}_{\text{epi}} = \text{interpolate\_basis}(\mathbf{Q}_{\text{LV,wall}}, \mathbf{Q}_{\text{RV,wall}}, \phi_{\text{BiV}}^{\text{interp}})$

7. **Interpolate from endocardium to epicardium**:
   - $\mathbf{Q} = \text{interpolate\_basis}(\mathbf{Q}_{\text{sep}}, \mathbf{Q}_{\text{epi}}, \phi_{\text{epi,trans}})$
   - Extract: $\mathbf{f} = \mathbf{Q}[:, 0]$, $\mathbf{n} = \mathbf{Q}[:, 1]$, $\mathbf{s} = \mathbf{Q}[:, 2]$




### Modified algorithm steps
In the original Doste paper, only one transmural field is mentioned ($\phi_{\text{Trans}}$) which is used to define the transmural gradient $\mathbf g_t$. Piersanti et al (2021) adds a second transmural field (see Fig. 4) ($\phi_{\text{BiV}}$) to differentiate between LV and RV. We use these two plus the individual transmural fields ($\phi_{\text{LV,trans}}$) and ($\phi_{\text{RV,trans}}$) which provides a smoother transition between LV and RV. 


# On the convention of the orthogonal basis and angles

Different papers use different conventions to define the orthogonal circumferential, longitudinal, and transmural basis (see Table below).

| Method | Longitudinal | Transmural | sign(alpha) endo/epi | sign(beta) endo/epi |
|--------|--------------|------------|----------------------|---------------------|
| Bayer | apex -> base | endo -> epi | +/- | -/+ |
| Doste | base -> apex | epi -> endo (rv) endo -> epi (lv) | -/+ (lv); +/- (rv) | -/+ |
| Piersanti (Doste) | apex -> base | endo -> epi (rv) epi -> endo (lv) | -/+ (lv); +/- (rv) in text section 2.6; +/- (lv); +/- (rv) in eq. 8 | -/+ |

For coherency, for all methods and for all chambers, we consider the transmural vector $\mathbf e_t$ pointing outwards (endo to epi), the longitudinal vector $\mathbf e_\ell$ pointint towards the base (apex to base), and the circumferential vector $\mathbf e_c=\mathbf e_\ell \times \mathbf e_t$. This way, when using literature angle values, we consider,

| Method | Longitudinal | Transmural | sign(alpha) endo/epi | sign(beta) endo/epi |
|--------|--------------|------------|----------------------|---------------------|
| FibGen | apex -> base | endo -> epi | +/- | -/+ |


# References
1. Bayer, J. D., Blake, R. C., Plank, G., & Trayanova, N. A. (2012). A Novel Rule-Based Algorithm for Assigning Myocardial Fiber Orientation to Computational Heart Models. Annals of Biomedical Engineering, 40(10), 2243–2254. https://doi.org/10.1007/s10439-012-0593-5
2. Doste, R., Soto‐Iglesias, D., Bernardino, G., Alcaine, A., Sebastian, R., Giffard‐Roisin, S., Sermesant, M., Berruezo, A., Sanchez‐Quintana, D., & Camara, O. (2019). A rule‐based method to model myocardial fiber orientation in cardiac biventricular geometries with outflow tracts. International Journal for Numerical Methods in Biomedical Engineering, 35(4). https://doi.org/10.1002/cnm.3185
3. Piersanti, R., Africa, P. C., Fedele, M., Vergara, C., Dedè, L., Corno, A. F., & Quarteroni, A. (2021). Modeling cardiac muscle fibers in ventricular and atrial electrophysiology simulations. Computer Methods in Applied Mechanics and Engineering, 373, 113468. https://doi.org/10.1016/j.cma.2020.113468