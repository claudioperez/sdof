pub struct State {
    displacement: f64,
    velocity: f64,
    acceleration: f64,
}

impl State {
    fn new() -> State {
        State {
            displacement: 0.0,
            velocity: 0.0,
            acceleration: 0.0,
        }
    }

    pub fn get(&self) -> (f64, f64, f64) {
        (self.displacement, self.velocity, self.acceleration)
    }

    pub fn set(&mut self, displacement: f64, velocity: f64, acceleration: f64) {
        self.displacement = displacement;
        self.velocity = velocity;
        self.acceleration = acceleration;
    }
}

pub trait Oscillator {
//  fn new(mass: f64, stiffness: f64, damping: f64) -> Box<dyn Oscillator>;
    fn get_state(&self) -> &State;
    fn set_state(&mut self, displacement: f64, velocity: f64, acceleration: f64);
    fn get_force(&self) -> f64;
    fn get_displ_tangent(&self) -> f64;
    fn get_veloc_tangent(&self) -> f64;
    fn get_accel_tangent(&self) -> f64;
}
pub struct DuffingOscillator {
    mass: f64,
    stiffness: f64,
    damping: f64,
    nonlinearity: f64,
    state: State,
}

impl DuffingOscillator {
    pub fn new(mass: f64, stiffness: f64, damping: f64, nonlinearity: f64) -> DuffingOscillator {
        DuffingOscillator {
            mass,
            stiffness,
            damping,
            nonlinearity,
            state: State::new(),
        }
    }
}

impl Oscillator for DuffingOscillator {
//  fn new(mass: f64, stiffness: f64, damping: f64) -> Box<dyn Oscillator> {
//      Box::new(DuffingOscillator {
//          mass,
//          stiffness,
//          damping,
//          nonlinearity: 1.0, // This is the nonlinearity parameter for the Duffing oscillator
//          state: State::new(),
//      })
//  }

    fn get_state(&self) -> &State {
        &self.state
    }
    
    fn set_state(&mut self, displacement: f64, velocity: f64, acceleration: f64) {
        self.state.set(displacement, velocity, acceleration);
    }

    fn get_force(&self) -> f64 {
        -self.stiffness * self.state.displacement - self.nonlinearity * self.state.displacement.powi(3)
    }

    fn get_displ_tangent(&self) -> f64 {
        -self.stiffness - 3.0 * self.nonlinearity * self.state.displacement.powi(2)
    }

    fn get_veloc_tangent(&self) -> f64 {
        -self.damping
    }

    fn get_accel_tangent(&self) -> f64 {
        self.mass
    }
}

pub struct LinearOscillator {
    mass: f64,
    stiffness: f64,
    damping: f64,
    state: State,
}

impl LinearOscillator {
    pub fn new(mass: f64, stiffness: f64, damping: f64) -> LinearOscillator {
//      Box::new(
        LinearOscillator {
            mass,
            stiffness,
            damping,
            state: State::new(),
        }
//  )
    }
}

impl Oscillator for LinearOscillator {

    fn get_state(&self) -> &State {
        &self.state
    }    
    fn set_state(&mut self, displacement: f64, velocity: f64, acceleration: f64) {
        self.state.set(displacement, velocity, acceleration);
    }

    fn get_force(&self) -> f64 {
        -self.stiffness * self.state.displacement
    }

    fn get_displ_tangent(&self) -> f64 {
        -self.stiffness
    }

    fn get_veloc_tangent(&self) -> f64 {
        -self.damping
    }

    fn get_accel_tangent(&self) -> f64 {
        self.mass
    }
}


pub struct VanDerPolOscillator {
    mass: f64,
    stiffness: f64,
    damping: f64,
    nonlinearity: f64,
    state: State,
}

impl VanDerPolOscillator {
    pub fn new(mass: f64, stiffness: f64, damping: f64, nonlinearity: f64) -> VanDerPolOscillator {
        VanDerPolOscillator {
            mass,
            stiffness,
            damping,
            nonlinearity,
            state: State::new(),
        }
    }
}

impl Oscillator for VanDerPolOscillator {
//  fn new(mass: f64, stiffness: f64, damping: f64) -> Box<dyn Oscillator> {
//      Box::new(VanDerPolOscillator {
//          mass,
//          stiffness,
//          damping,
//          nonlinearity: 1.0, // This is the nonlinearity parameter for the Van der Pol oscillator
//          state: State::new(),
//      })
//  }

    fn get_state(&self) -> &State {
        &self.state
    }
    
    fn set_state(&mut self, displacement: f64, velocity: f64, acceleration: f64) {
        self.state.set(displacement, velocity, acceleration);
    }

    fn get_force(&self) -> f64 {
        self.nonlinearity * (1.0 - self.state.displacement.powi(2)) * self.state.velocity - self.state.displacement
    }

    fn get_displ_tangent(&self) -> f64 {
        -2.0 * self.nonlinearity * self.state.displacement * self.state.velocity - 1.0
    }

    fn get_veloc_tangent(&self) -> f64 {
        self.nonlinearity * (1.0 - self.state.displacement.powi(2))
    }

    fn get_accel_tangent(&self) -> f64 {
        self.mass
    }
}

