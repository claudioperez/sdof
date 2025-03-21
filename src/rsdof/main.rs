mod oscillator;
mod integrator;

use std::io::{self, BufRead};
use std::env;

use oscillator::Oscillator;
use integrator::Integrator;
use integrator::NewmarkIntegrator;
use oscillator::LinearOscillator;
use oscillator::DuffingOscillator;
use oscillator::VanDerPolOscillator;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut oscillator: Box<dyn Oscillator>;

    match args.len() {
        3 => {
            let mass: f64 = args[1].parse().unwrap();
            let stiffness: f64 = args[2].parse().unwrap();
            oscillator = Box::new(LinearOscillator::new(mass, 0.0, stiffness));
        },
        4 => {
            let mass: f64 = args[1].parse().unwrap();
            let damping: f64 = args[2].parse().unwrap();
            let stiffness: f64 = args[3].parse().unwrap();
            oscillator = Box::new(LinearOscillator::new(mass, damping, stiffness));
        },
        _ => {
            let mass: f64 = args[1].parse().unwrap();
            let damping: f64 = args[2].parse().unwrap();
            let stiffness: f64 = args[3].parse().unwrap();

            match args[4].as_str() {
                "--van" => {
                    let nonlinearity: f64 = args[5].parse().unwrap();
                    oscillator = Box::new(VanDerPolOscillator::new(mass, damping, stiffness, nonlinearity));
                },
                "--duf" => {
                    let nonlinearity: f64 = args[5].parse().unwrap();
                    oscillator = Box::new(DuffingOscillator::new(mass, damping, stiffness, nonlinearity));
                },
                _ => {
                    eprintln!("Invalid argument: {}", args[4]);
                    return;
                },
            }
        },
    }

    if args.len() > 4 {
    }

//  let integrator: Box<dyn Integrator> = Box::new(NewmarkIntegrator::new(0.5, 0.25, 0.1));
    let integrator = NewmarkIntegrator::new(0.5, 0.25, 0.1);
    let stdin = io::stdin();

    for line in stdin.lock().lines() {
        let force: f64 = line.unwrap().parse().unwrap();
        integrator.integrate(&mut *oscillator, force, 0.01);
        println!("{}", oscillator.get_state().get().0);
    }
}

