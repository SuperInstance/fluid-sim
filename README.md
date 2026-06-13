# 2D Fluid Simulation

**An Eulerian grid-based fluid simulator in Rust** implementing the Stam (1999) stable fluids algorithm with semi-Lagrangian advection, Gauss-Seidel diffusion, and pressure projection on a 32×32 grid.

## Why It Matters

Fluid simulation is fundamental to visual effects (smoke, fire, water), scientific computing (weather models, aerodynamics), and game physics. The "Stable Fluids" method published by Jos Stam in 1999 revolutionized the field by introducing an unconditionally stable advection scheme, enabling real-time simulation without tiny time steps. This implementation is a self-contained, dependency-free reference that teaches the core algorithm: density advection, viscous diffusion, and Helmholtz-Hodge decomposition via pressure projection.

## How It Works

The simulation runs on a 34×34 grid (32 interior cells + 1-cell boundary on each side) using a fixed `N=32` resolution. Each timestep executes four passes:

1. **Add source**: Density and velocity sources are injected via `x[i] += dt * s[i]` — **O(N²)**.

2. **Diffuse**: Solves the heat equation `x - x₀ = a·∇²x` using 20 iterations of Gauss-Seidel relaxation, where `a = dt·diff·N²`. Each iteration touches the 4-neighborhood of every cell — **O(N²)** per iteration.

3. **Advect**: Uses semi-Lagrangian backtracing. For each cell `(i,j)`, it traces backward through the velocity field to find the source position, clamps to grid bounds, and bilinearly interpolates the previous density field. This is the key insight from Stam: tracing backward along characteristics guarantees no instability — **O(N²)**.

4. **Project**: Enforces incompressibility (∇·**u** = 0) by solving a Poisson equation for pressure. The divergence is computed via central differences, then 20 Gauss-Seidel iterations solve `∇²p = ∇·**u**/h`, and the pressure gradient is subtracted from the velocity field — **O(N²)** per iteration.

Boundary conditions (`set_bnd`) handle walls: for velocity components, the normal component reflects (negated), while tangential components copy. Density uses simple copying at boundaries.

## Quick Start

```bash
cargo run
```

This runs 20 simulation steps with a density/velocity source injected in the center grid for the first 5 steps, then watches it advect and diffuse:

```
Step  0: max density = 10.0100
Step  1: max density = 20.0200
...
Step 19: max density = 0.0003
Fluid sim complete.
```

## API

The library is a binary crate (`main.rs`). Key internal functions:

| Function | Description |
|---|---|
| `add_source(x, s, dt)` | Add source field to state field — **O(N²)** |
| `set_bnd(b, x)` | Apply boundary conditions (reflect/copy) |
| `diffuse(b, x, x0, diff, dt)` | Gauss-Seidel diffusion solver — **O(N²·iters)** |
| `advect(b, d, d0, u, v, dt)` | Semi-Lagrangian advection — **O(N²)** |
| `project(u, v, p, div)` | Pressure projection for incompressibility — **O(N²·iters)** |

## Architecture Notes

Part of the SuperInstance scientific computing collection. This is a standalone reference implementation. See the [Architecture Guide](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## License

MIT
