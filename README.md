# fluid-sim: 2D Eulerian Fluid Simulation

A from-scratch implementation of the **Stam (1999) Stable Fluids** algorithm — a grid-based (Eulerian) fluid solver that is unconditionally stable via a semi-Lagrangian advection scheme. The simulation evolves a velocity field **u** = (u, v) and a density field **ρ** through the Navier–Stokes equations for incompressible flow.

## Why It Matters

Fluid simulation is the backbone of visual effects, weather modeling, and aerodynamics. The Stable Fluids method was revolutionary because it decoupled simulation stability from timestep size — previous explicit schemes required `dt` small enough to satisfy CFL conditions that made interactive simulation impossible. This implementation demonstrates every core operator in a computational fluid dynamics (CFD) solver in under 200 lines of readable Rust.

## How It Works

The solver advances the incompressible Navier–Stokes equations:

```
∂u/∂t = -(u·∇)u + ν∇²u - ∇p + f
∇·u = 0
```

Each timestep applies four operators sequentially:

### 1. Diffusion (Implicit Gauss–Seidel)

Viscous diffusion solves `x' = x + dt·ν·∇²x` via 20 iterations of Gauss–Seidel relaxation on the discretized Laplacian:

```
x[i,j] = (x₀[i,j] + a·(x[i-1,j] + x[i+1,j] + x[i,j-1] + x[i,j+1])) / (1 + 4a)
```

where `a = dt·ν·N²`. This is **O(N²·k)** per call where k = 20 iterations.

### 2. Advection (Semi-Lagrangian)

Instead of tracking particles, the semi-Lagrangian scheme back-traces from each grid cell through the velocity field and interpolates:

```
x_new[i,j] = bilinear(d₀, x_back, y_back)
```

This is **O(N²)** and unconditionally stable — the key insight from Stam (1999).

### 3. Projection (Pressure Solve)

To enforce incompressibility (∇·u = 0), we solve a Poisson equation for pressure:

```
∇²p = ∇·u / h
```

via 20 Gauss–Seidel iterations, then subtract the pressure gradient from velocity. This is **O(N²·k)** with k = 20.

### 4. Boundary Conditions

`set_bnd` enforces reflective (no-slip for tangential) or free-slip boundary conditions on the grid edges depending on the field type (velocity component vs. density).

### Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Diffusion | O(N²·k) | O(N²) |
| Advection | O(N²) | O(N²) |
| Projection | O(N²·k) | O(N²) |
| Full step | O(N²·k) | O(N²) |

With N = 32 and k = 20: ~40K cell-updates per operator per step.

## Quick Start

```bash
cargo run
```

Outputs density maxima over 20 timesteps:

```
Step  0: max density = 10.0000
Step  1: max density = 19.8000
...
Step 19: max density = 3.2145
```

## API

The simulation is contained in `src/main.rs` with these core functions:

| Function | Signature | Purpose |
|----------|-----------|---------|
| `add_source` | `(&mut [f64; SIZE], &[f64; SIZE], f64)` | Add source terms scaled by dt |
| `diffuse` | `(i32, &mut [f64; SIZE], &[f64; SIZE], f64, f64)` | Implicit diffusion via Gauss–Seidel |
| `advect` | `(i32, &mut [f64; SIZE], &[f64; SIZE], &[f64; SIZE], &[f64; SIZE], f64)` | Semi-Lagrangian advection |
| `project` | `(&mut [f64; SIZE], &mut [f64; SIZE], &mut [f64; SIZE], &mut [f64; SIZE])` | Pressure projection for incompressibility |
| `set_bnd` | `(i32, &mut [f64; SIZE])` | Apply boundary conditions |

Grid resolution: `N = 32` (internal cells), total array size `(N+2)² = 1156`.

## Architecture Notes

This implementation fits into the **γ + η = C** framework as a concrete **γ (gamma)** module — a physics solver that produces deterministic, testable outputs. It can serve as the computational core for an agent that reasons about physical systems. The density field is the observable state; the velocity field is the hidden dynamic that the agent must infer or control.

### Numerical Parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `N` | 32 | Grid resolution (32×32 internal cells) |
| `dt` | 0.1 | Timestep (seconds) |
| `diff` | 0.0001 | Density diffusion coefficient |
| `visc` | 0.0 | Velocity viscosity (inviscid) |
| Gauss–Seidel iterations | 20 | Per solve for diffusion and projection |

With N = 32 and k = 20: ~40K cell-updates per operator per step.

## References

- Stam, J. (1999). *Stable Fluids*. Proceedings of SIGGRAPH '99, pp. 121–128.
- Bridson, R. (2015). *Fluid Simulation for Computer Graphics* (2nd ed.). CRC Press.
- Fedkiw, R., Stam, J., & Jensen, H. W. (2001). *Visual Simulation of Smoke*. ACM TOG 20(3).

## License

MIT
