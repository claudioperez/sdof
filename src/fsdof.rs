#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#![register_tool(c2rust)]
#![feature(register_tool)]
extern "C" {
    fn abs(_: f64) -> f64;
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct generalized_alpha {
    pub alpha_m: f64,
    pub alpha_f: f64,
    pub beta: f64,
    pub gamma: f64,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct SDOF_Peaks {
    pub max_displ: f64,
    pub max_accel: f64,
    pub time_max_accel: f64,
}
#[no_mangle]
pub static mut CONF: generalized_alpha = {
    let mut init = generalized_alpha {
        alpha_m: 1.0f64,
        alpha_f: 1.0f64,
        beta: 0.25f64,
        gamma: 0.5f64,
    };
    init
};

pub unsafe extern "C" fn fsdof_peaks(
    mut conf: *mut generalized_alpha,
    mut M: f64,
    mut C: f64,
    mut K: f64,
    mut scale: f64,
    mut n: i32,
    mut p: *mut f64,
    mut dt: f64,
    mut response: *mut SDOF_Peaks,
) -> i32 {
    conf = &mut CONF;
    let gamma: f64 = (*conf).gamma;
    let beta: f64 = (*conf).beta;
    let alpha_m: f64 = (*conf).alpha_m;
    let alpha_f: f64 = (*conf).alpha_f;
    let c1: f64 = 1.0f64;
    let c2: f64 = gamma / (beta * dt);
    let c3: f64 = 1.0f64 / (beta * dt * dt);
    let a1: f64 = 1.0f64 - gamma / beta;
    let a2: f64 = dt * (1.0f64 - 0.5f64 * gamma / beta);
    let a3: f64 = -1.0f64 / (beta * dt);
    let a4: f64 = 1.0f64 - 0.5f64 / beta;
    let ki: f64 = alpha_f * c1 * K + alpha_f * c2 * C + alpha_m * c3 * M;
    let mut time: f64 = 0.0f64;
    let mut ua: f64 = 0.;
    let mut va: f64 = 0.;
    let mut aa: f64 = 0.;
    let mut u: [f64; 2] = [0.; 2];
    let mut v: [f64; 2] = [0.; 2];
    let mut a: [f64; 2] = [0.; 2];
    let mut i: i32 = 0 as i32;
    let mut past: i32 = 1 as i32;
    let mut pres: i32 = 0 as i32;
    u[pres as usize] = 0.0f64;
    v[pres as usize] = 0.0f64;
    a[pres as usize] = (*p.offset(i as isize) - C*v[pres as usize] - K*u[pres as usize])/M;
    i = 1 as i32;
    while i < n {
        past = (past == 0) as i32;
        pres = (pres == 0) as i32;
        u[pres as usize] = u[past as usize];
        v[pres as usize] = a1 * v[past as usize] + a2 * a[past as usize];
        a[pres as usize] = a4 * a[past as usize] + a3 * v[past as usize];
        va = (1 as i32 as f64 - alpha_f) * v[past as usize]
            + alpha_f * v[pres as usize];
        aa = (1 as i32 as f64 - alpha_m) * a[past as usize]
            + alpha_m * a[pres as usize];
        let mut pi: f64 = scale * *p.offset(i as isize) - C * va - M * aa
            - K * u[pres as usize];
        let mut du: f64 = pi / ki;
        u[pres as usize] += du;
        v[pres as usize] += c2 * du;
        a[pres as usize] += c3 * du;
        if f64::abs(u[pres as usize]) > (*response).max_displ {
            (*response).max_displ = f64::abs(u[pres as usize]);
        }
        if f64::abs(a[pres as usize]) > (*response).max_accel {
            (*response).max_accel = f64::abs(a[pres as usize]);
            (*response).time_max_accel = i as f64 * dt;
        }
        i += 1;
    }
    return 1;
}

pub unsafe extern "C" fn fsdof_integrate(
    mut conf: *mut generalized_alpha,
    M: f64, C: f64, K: f64,
    mut scale: f64,
    mut n: i32,
    mut p: *mut f64,
    mut dt: f64,
    mut response: *mut f64,
) -> i32 {
    conf = &mut CONF;
    let gamma: f64 = (*conf).gamma;
    let beta: f64 = (*conf).beta;
    let alpha_m: f64 = (*conf).alpha_m;
    let alpha_f: f64 = (*conf).alpha_f;
    let c1: f64 = 1.0f64;
    let c2: f64 = gamma / (beta * dt);
    let c3: f64 = 1.0f64 / (beta * dt * dt);
    let a1: f64 = 1.0f64 - gamma / beta;
    let a2: f64 = dt * (1.0f64 - 0.5f64 * gamma / beta);
    let a3: f64 = -1.0f64 / (beta * dt);
    let a4: f64 = 1.0f64 - 0.5f64 / beta;
    let ki: f64 = alpha_f * c1 * K + alpha_f * c2 * C + alpha_m * c3 * M;
    let mut time: f64 = 0.0f64;
    let mut ua: f64 = 0.;
    let mut va: f64 = 0.;
    let mut aa: f64 = 0.;
    let mut u: *mut f64 = &mut *response.offset(0 as i32 as isize)
        as *mut f64;
    let mut v: *mut f64 = &mut *response.offset(n as isize)
        as *mut f64;
    let mut a: *mut f64 = &mut *response
        .offset((2 as i32 * n) as isize) as *mut f64;
    let mut i: i32 = 0 as i32;
    let past: i32 = -(1 as i32);
    let pres: i32 = 0 as i32;
    *a
        .offset(
            pres as isize,
        ) = (*p.offset(i as isize) - C * *v.offset(pres as isize)
        - K * *u.offset(pres as isize)) / M;
    i = 1 as i32;
    while i < n {
        u = u.offset(1);
        v = v.offset(1);
        a = a.offset(1);
        *u.offset(pres as isize) = *u.offset(past as isize);
        *v
            .offset(
                pres as isize,
            ) = a1 * *v.offset(past as isize) + a2 * *a.offset(past as isize);
        *a
            .offset(
                pres as isize,
            ) = a4 * *a.offset(past as isize) + a3 * *v.offset(past as isize);
        va = (1.0f64 - alpha_f) * *v.offset(past as isize)
            + alpha_f * *v.offset(pres as isize);
        aa = (1.0f64 - alpha_m) * *a.offset(past as isize)
            + alpha_m * *a.offset(pres as isize);
        let mut pi: f64 = scale * *p.offset(i as isize) - C * va - M * aa
            - K * *u.offset(pres as isize);
        let mut du: f64 = pi / ki;
        *u.offset(pres as isize) += du;
        *v.offset(pres as isize) += c2 * du;
        *a.offset(pres as isize) += c3 * du;
        i += 1;
    }
    return 1;
}

pub unsafe extern "C" fn fsdof_integrate2(
    mut conf: *mut generalized_alpha,
    M: f64, C: f64, K: f64,
    mut scale: f64,
    mut n: i32, mut p: *mut f64,
    mut dt: f64,
    mut response: *mut f64,
) -> i32 {
    conf = &mut CONF;
    let gamma: f64 = (*conf).gamma;
    let beta: f64 = (*conf).beta;
    let alpha_m: f64 = (*conf).alpha_m;
    let alpha_f: f64 = (*conf).alpha_f;
    let c1: f64 = 1.0f64;
    let c2: f64 = gamma / (beta * dt);
    let c3: f64 = 1.0f64 / (beta * dt * dt);
    let a1: f64 = 1.0f64 - gamma / beta;
    let a2: f64 = dt * (1.0f64 - 0.5f64 * gamma / beta);
    let a3: f64 = -1.0f64 / (beta * dt);
    let a4: f64 = 1.0f64 - 0.5f64 / beta;
    let ki: f64 = alpha_f * c1 * K + alpha_f * c2 * C + alpha_m * c3 * M;
    let mut time: f64 = 0.0f64;
    let mut ua: f64 = 0.;
    let mut va: f64 = 0.;
    let mut aa: f64 = 0.;
    let mut u: *mut f64 = &mut *response.offset(0 as i32 as isize)
        as *mut f64;
    let mut v: *mut f64 = &mut *response.offset(1 as i32 as isize)
        as *mut f64;
    let mut a: *mut f64 = &mut *response.offset(2 as i32 as isize)
        as *mut f64;
    let mut i: i32 = 0 as i32;
    let past: i32 = -(3 as i32);
    let pres: i32 = 0 as i32;
    *a.offset(
            pres as isize,
        ) = (*p.offset(i as isize) - C * *v.offset(pres as isize)
        - K * *u.offset(pres as isize)) / M;

    i = 1 as i32;
    while i < n {
        u = u.offset(3 as i32 as isize);
        v = v.offset(3 as i32 as isize);
        a = a.offset(3 as i32 as isize);
        *u.offset(pres as isize) = *u.offset(past as isize);
        *v.offset(
                pres as isize,
            ) = a1 * *v.offset(past as isize) + a2 * *a.offset(past as isize);
        *a.offset(
                pres as isize,
            ) = a4 * *a.offset(past as isize) + a3 * *v.offset(past as isize);
        va = (1.0f64 - alpha_f) * *v.offset(past as isize)
            + alpha_f * *v.offset(pres as isize);
        aa = (1.0f64 - alpha_m) * *a.offset(past as isize)
            + alpha_m * *a.offset(pres as isize);
        let mut pi: f64 = scale * *p.offset(i as isize) - C * va - M * aa
            - K * *u.offset(pres as isize);
        let mut du: f64 = pi / ki;
        *u.offset(pres as isize) += du;
        *v.offset(pres as isize) += c2 * du;
        *a.offset(pres as isize) += c3 * du;
        i += 1;
    }
    return 1;
}
