use crate::oscillator::Oscillator;


pub trait Integrator {
    fn integrate(&self, oscillator: &mut dyn Oscillator, force: f64, dt: f64);
}

pub struct NewmarkIntegrator {
    gamma: f64,
    beta: f64,
    alpha: f64,
}

impl Integrator for NewmarkIntegrator {
    fn new(gamma: f64, beta: f64, alpha: f64) -> NewmarkIntegrator {
        NewmarkIntegrator { gamma, beta, alpha }
    }

    fn integrate(&self, oscillator: &mut dyn Oscillator, force: f64, dt: f64) {
        let (mut u, mut v, mut a) = oscillator.get_state().get();

        // Newton-Raphson method for handling nonlinearity
        let mut delta_u;
        let mut residual;
        let mut jacobian;
        let tolerance = 1e-6;
        let max_iterations = 100;

        for _ in 0..max_iterations {
            let response = oscillator.get_force();
            let displ_tangent = oscillator.get_displ_tangent();
            let veloc_tangent = oscillator.get_veloc_tangent();
            let accel_tangent = oscillator.get_accel_tangent();

            let effective_stiffness = displ_tangent - self.gamma / (self.beta * dt) * veloc_tangent + accel_tangent / (self.beta * dt * dt);
            let effective_force = force - response + veloc_tangent * (self.gamma / (self.beta * dt) * u - v) + accel_tangent * (1.0 / (self.beta * dt) * u + v / self.beta);

            residual = effective_force - accel_tangent * a;
            jacobian = effective_stiffness;

            delta_u = -residual / jacobian;

            if delta_u.abs() < tolerance {
                break;
            }

            u += delta_u;
        }

        let response      = oscillator.get_force();
        let veloc_tangent = oscillator.get_veloc_tangent();
        let accel_tangent = oscillator.get_accel_tangent();

        let effective_force = force - response + veloc_tangent * (self.gamma / (self.beta * dt) * u - v) + accel_tangent * (1.0 / (self.beta * dt) * u + v / self.beta);

        a = effective_force / accel_tangent;
        v += dt * ((1.0 - self.gamma) * a + self.gamma * a);
        u += dt * v + dt * dt * (0.5 - self.beta) * a;

        oscillator.set_state(u, v, a);
    }
}

