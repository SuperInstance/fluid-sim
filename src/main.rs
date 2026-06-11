/// Simplified 2D fluid simulation using a grid-based approach (Eulerian)
/// Implements basic advection and diffusion on a density field

const N: usize = 32;
const SIZE: usize = (N + 2) * (N + 2);

fn ix(i: usize, j: usize) -> usize { i + (N + 2) * j }

fn add_source(x: &mut [f64; SIZE], s: &[f64; SIZE], dt: f64) {
    for i in 0..SIZE { x[i] += dt * s[i]; }
}

fn set_bnd(b: i32, x: &mut [f64; SIZE]) {
    for i in 1..=N {
        x[ix(0, i)]     = if b == 1 { -x[ix(1, i)] } else { x[ix(1, i)] };
        x[ix(N+1, i)]   = if b == 1 { -x[ix(N, i)] } else { x[ix(N, i)] };
        x[ix(i, 0)]     = if b == 2 { -x[ix(i, 1)] } else { x[ix(i, 1)] };
        x[ix(i, N+1)]   = if b == 2 { -x[ix(i, N)] } else { x[ix(i, N)] };
    }
    x[ix(0, 0)]       = 0.5 * (x[ix(1, 0)] + x[ix(0, 1)]);
    x[ix(0, N+1)]     = 0.5 * (x[ix(1, N+1)] + x[ix(0, N)]);
    x[ix(N+1, 0)]     = 0.5 * (x[ix(N, 0)] + x[ix(N+1, 1)]);
    x[ix(N+1, N+1)]   = 0.5 * (x[ix(N, N+1)] + x[ix(N+1, N)]);
}

fn diffuse(b: i32, x: &mut [f64; SIZE], x0: &[f64; SIZE], diff: f64, dt: f64) {
    let a = dt * diff * (N as f64) * (N as f64);
    for _ in 0..20 {
        for j in 1..=N {
            for i in 1..=N {
                x[ix(i,j)] = (x0[ix(i,j)] + a * (
                    x[ix(i-1,j)] + x[ix(i+1,j)] +
                    x[ix(i,j-1)] + x[ix(i,j+1)]
                )) / (1.0 + 4.0 * a);
            }
        }
        set_bnd(b, x);
    }
}

fn advect(b: i32, d: &mut [f64; SIZE], d0: &[f64; SIZE], u: &[f64; SIZE], v: &[f64; SIZE], dt: f64) {
    let dt0 = dt * N as f64;
    for j in 1..=N {
        for i in 1..=N {
            let mut x = i as f64 - dt0 * u[ix(i,j)];
            let mut y = j as f64 - dt0 * v[ix(i,j)];
            x = x.max(0.5).min(N as f64 + 0.5);
            y = y.max(0.5).min(N as f64 + 0.5);
            let i0 = x as usize; let i1 = i0 + 1;
            let j0 = y as usize; let j1 = j0 + 1;
            let s1 = x - i0 as f64; let s0 = 1.0 - s1;
            let t1 = y - j0 as f64; let t0 = 1.0 - t1;
            d[ix(i,j)] = s0*(t0*d0[ix(i0,j0)] + t1*d0[ix(i0,j1)])
                        + s1*(t0*d0[ix(i1,j0)] + t1*d0[ix(i1,j1)]);
        }
    }
    set_bnd(b, d);
}

fn project(u: &mut [f64; SIZE], v: &mut [f64; SIZE], p: &mut [f64; SIZE], div: &mut [f64; SIZE]) {
    let h = 1.0 / N as f64;
    for j in 1..=N {
        for i in 1..=N {
            div[ix(i,j)] = -0.5 * h * (
                u[ix(i+1,j)] - u[ix(i-1,j)] + v[ix(i,j+1)] - v[ix(i,j-1)]
            );
            p[ix(i,j)] = 0.0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    for _ in 0..20 {
        for j in 1..=N {
            for i in 1..=N {
                p[ix(i,j)] = (div[ix(i,j)] + p[ix(i-1,j)] + p[ix(i+1,j)]
                    + p[ix(i,j-1)] + p[ix(i,j+1)]) / 4.0;
            }
        }
        set_bnd(0, p);
    }
    for j in 1..=N {
        for i in 1..=N {
            u[ix(i,j)] -= 0.5 * (p[ix(i+1,j)] - p[ix(i-1,j)]) / h;
            v[ix(i,j)] -= 0.5 * (p[ix(i,j+1)] - p[ix(i,j-1)]) / h;
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
}

fn main() {
    let dt = 0.1;
    let diff = 0.0001;
    let visc = 0.0;

    let mut dens = [0.0f64; SIZE]; let mut dens_prev = [0.0f64; SIZE];
    let mut u = [0.0f64; SIZE];    let mut u_prev = [0.0f64; SIZE];
    let mut v = [0.0f64; SIZE];    let mut v_prev = [0.0f64; SIZE];

    for step in 0..20 {
        // Add source in center
        if step < 5 {
            for j in 14..=18 {
                for i in 14..=18 {
                    dens_prev[ix(i,j)] = 100.0;
                    u_prev[ix(i,j)] = 5.0;
                    v_prev[ix(i,j)] = 2.0;
                }
            }
        }

        // Velocity step
        add_source(&mut u, &u_prev, dt);
        add_source(&mut v, &v_prev, dt);

        std::mem::swap(&mut u_prev, &mut u);
        diffuse(1, &mut u, &u_prev, visc, dt);
        std::mem::swap(&mut v_prev, &mut v);
        diffuse(2, &mut v, &v_prev, visc, dt);
        project(&mut u, &mut v, &mut u_prev, &mut v_prev);

        std::mem::swap(&mut u_prev, &mut u);
        std::mem::swap(&mut v_prev, &mut v);
        {
            let u_src = u_prev;
            let v_src = v_prev;
            advect(1, &mut u, &mut u_prev, &u_src, &v_src, dt);
            advect(2, &mut v, &mut v_prev, &u_src, &v_src, dt);
        }
        project(&mut u, &mut v, &mut u_prev, &mut v_prev);

        // Density step
        add_source(&mut dens, &dens_prev, dt);
        std::mem::swap(&mut dens_prev, &mut dens);
        diffuse(0, &mut dens, &mut dens_prev, diff, dt);
        std::mem::swap(&mut dens_prev, &mut dens);
        advect(0, &mut dens, &mut dens_prev, &mut u, &mut v, dt);

        dens_prev = [0.0; SIZE];
        u_prev = [0.0; SIZE];
        v_prev = [0.0; SIZE];

        // Report max density
        let max_d = dens.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("Step {step:2}: max density = {max_d:.4}");
    }
    println!("Fluid sim complete.");
}
