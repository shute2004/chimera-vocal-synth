//! Diagnostic only: excite a linearized tract from either endpoint and measure
//! pressure at the other. This is not reverse speech synthesis or inversion.
use chimera::airway::network::LinearElements;
use chimera::anatomy::human::{ArticKey, BodyMix};
use chimera::phonetics::planner::vowel_artic;
use chimera::physics::{C_SOUND, RHO};
use chimera::voice::human_voice::HumanVoice;
use std::ops::{Add, Div, Mul};

#[derive(Copy, Clone, Default, Debug)]
struct Z {
    r: f64,
    i: f64,
}
impl Z {
    fn new(r: f64, i: f64) -> Self {
        Self { r, i }
    }
    fn abs(self) -> f64 {
        self.r.hypot(self.i)
    }
}
impl Add for Z {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Z::new(self.r + o.r, self.i + o.i)
    }
}
impl Mul for Z {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Z::new(self.r * o.r - self.i * o.i, self.r * o.i + self.i * o.r)
    }
}
impl Div for Z {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let d = o.r * o.r + o.i * o.i;
        Z::new(
            (self.r * o.r + self.i * o.i) / d,
            (self.i * o.r - self.r * o.i) / d,
        )
    }
}
fn response(el: &LinearElements, inject: usize, observe: usize, cut: usize, hz: f64) -> f64 {
    let n = el.n_nodes;
    let w = std::f64::consts::TAU * hz;
    let mut adj = vec![Vec::<(usize, usize)>::new(); n];
    for j in 0..el.n_branches {
        if !el.open[j] || j == cut {
            continue;
        }
        adj[el.a[j]].push((el.b[j], j));
        adj[el.b[j]].push((el.a[j], j));
    }
    let mut parent = vec![usize::MAX; n];
    let mut edge = vec![usize::MAX; n];
    let mut order = vec![inject];
    parent[inject] = inject;
    let mut pos = 0;
    while pos < order.len() {
        let k = order[pos];
        pos += 1;
        for &(other, j) in &adj[k] {
            if parent[other] != usize::MAX {
                continue;
            }
            parent[other] = k;
            edge[other] = j;
            order.push(other);
        }
    }
    assert!(parent[observe] != usize::MAX, "Endpoints disconnected");
    let mut y = vec![Z::default(); n];
    for &k in &order {
        y[k] = Z::new(if el.fixed[k] { 1e30 } else { 0.0 }, w * el.c[k]);
        if el.wall_area[k] > 0.0 {
            let z = Z::new(el.wall_r[k], w * el.wall_m[k] - el.wall_k[k] / w);
            y[k] = y[k] + Z::new(el.wall_area[k], 0.0) / z;
        }
    }
    for (r, &node) in el.rad_nodes.iter().enumerate() {
        let a = el.rad_area[r];
        let lr = (8.0 * RHO) / (3.0 * std::f64::consts::PI * (std::f64::consts::PI * a).sqrt());
        let rr = (128.0 * RHO * C_SOUND) / (9.0 * std::f64::consts::PI.powi(2) * a);
        y[node] = y[node] + Z::new(1.0 / rr, -1.0 / (w * lr));
    }
    let mut denom = vec![Z::default(); n];
    for &k in order.iter().skip(1).rev() {
        let j = edge[k];
        let z = Z::new(el.r[j], w * el.l[j]);
        denom[k] = Z::new(1.0, 0.0) + z * y[k];
        let eq = y[k] / denom[k];
        y[parent[k]] = y[parent[k]] + eq;
    }
    let mut p = vec![Z::default(); n];
    p[inject] = Z::new(1.0, 0.0) / y[inject];
    for &k in order.iter().skip(1) {
        p[k] = p[parent[k]] / denom[k];
    }
    p[observe].abs()
}
#[test]
fn lip_to_glottis_matches_glottis_to_lip_for_fixed_tract() {
    let mut signatures = Vec::new();
    for ipa in ["a", "i"] {
        let body = BodyMix::new(0.45, 0.0);
        let mut v = HumanVoice::new(48000.0, body.body(), 2);
        v.set_articulation(&vowel_artic(ipa, body).unwrap(), None);
        v.set_gate(1.0);
        v.render(&mut vec![0.0f32; 9600]);
        let el = v.airway.net.linear_elements();
        let lip = v.airway.oral.last_node();
        let glottis = v.airway.pharynx.first_node();
        let cut = v.airway.glottis_port;
        let mut levels = Vec::new();
        for hz in [300.0, 600.0, 900.0, 1200.0] {
            let forward = response(&el, glottis, lip, cut, hz);
            let reverse = response(&el, lip, glottis, cut, hz);
            let rel = (forward - reverse).abs() / forward.max(reverse).max(1e-30);
            assert!(rel < 1e-8, "{ipa} {hz} reciprocity error {rel}");
            levels.push(format!("{:.5e}", forward));
        }
        println!(
            "{ipa}: |P(lips|glottis)| = |P(glottis|lips)| at [300,600,900,1200]Hz: {}",
            levels.join(", ")
        );
        signatures.push(levels);
    }
    assert_ne!(
        signatures[0], signatures[1],
        "Different tract shapes should have different responses"
    );
}

/// Probe the lips with known volume-velocity amplitudes and invert the four
/// observed glottal pressures for one articulator. A known vowel is assumed.
#[test]
fn recover_jaw_from_external_lip_probe() {
    fn probe(jaw_shift: f64) -> [f64; 12] {
        let body = BodyMix::new(0.45, 0.0);
        let mut v = HumanVoice::new(48000.0, body.body(), 2);
        let mut artic = vowel_artic("a", body).unwrap();
        let jaw = artic.get(ArticKey::Jaw).unwrap();
        artic.set(ArticKey::Jaw, jaw + jaw_shift);
        v.set_articulation(&artic, None);
        v.set_gate(1.0);
        v.render(&mut vec![0.0f32; 9600]);
        let el = v.airway.net.linear_elements();
        [
            250.0, 400.0, 600.0, 800.0, 1100.0, 1400.0, 1700.0, 2000.0, 2400.0, 2800.0, 3300.0,
            3800.0,
        ]
        .map(|hz| {
            response(
                &el,
                v.airway.oral.last_node(),
                v.airway.pharynx.first_node(),
                v.airway.glottis_port,
                hz,
            )
            .ln()
        })
    }
    let target_jaw = 0.27;
    let target = probe(target_jaw);
    let error = |values: [f64; 12]| -> f64 {
        (0..12)
            .map(|i| (values[i] - target[i]).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let mut jaw = -0.30;
    let before = error(probe(jaw));
    let mut score = before;
    for _ in 0..16 {
        let current = probe(jaw);
        let plus = probe((jaw + 0.015).min(0.7));
        let minus = probe((jaw - 0.015).max(-0.7));
        let deriv = std::array::from_fn::<f64, 12, _>(|i| {
            (plus[i] - minus[i]) / ((jaw + 0.015).min(0.7) - (jaw - 0.015).max(-0.7))
        });
        let grad = (0..12)
            .map(|i| (current[i] - target[i]) * deriv[i])
            .sum::<f64>();
        let curvature = (0..12).map(|i| deriv[i] * deriv[i]).sum::<f64>() + 0.01;
        let direction = (-grad / curvature).clamp(-0.2, 0.2);
        let mut improved = false;
        for alpha in [1.0, 0.5, 0.25] {
            let guess = (jaw + alpha * direction).clamp(-0.7, 0.7);
            let loss = error(probe(guess));
            if loss < score - 1e-8 {
                jaw = guess;
                score = loss;
                improved = true;
                break;
            }
        }
        if !improved {
            break;
        }
    }
    println!(
        "mouth-probe inverse: true jaw shift={target_jaw:.5}, recovered={jaw:.5}, loss {before:.6} -> {score:.6}"
    );
    assert!(
        score < before * 0.01,
        "Inverse did not reduce response mismatch"
    );
    assert!(
        (jaw - target_jaw).abs() < 0.03,
        "Recovered jaw differs from hidden setting"
    );
}

#[test]
fn recover_jaw_and_tongue_height_from_lip_probes() {
    fn probe(jaw_shift: f64, tongue_shift: f64) -> [f64; 12] {
        let body = BodyMix::new(0.45, 0.0);
        let mut v = HumanVoice::new(48000.0, body.body(), 2);
        let mut artic = vowel_artic("a", body).unwrap();
        let jaw = artic.get(ArticKey::Jaw).unwrap();
        let tongue = artic.get(ArticKey::TbY).unwrap();
        artic.set(ArticKey::Jaw, jaw + jaw_shift);
        artic.set(ArticKey::TbY, tongue + tongue_shift);
        v.set_articulation(&artic, None);
        v.set_gate(1.0);
        v.render(&mut vec![0.0f32; 9600]);
        let el = v.airway.net.linear_elements();
        [
            250.0, 400.0, 600.0, 800.0, 1100.0, 1400.0, 1700.0, 2000.0, 2400.0, 2800.0, 3300.0,
            3800.0,
        ]
        .map(|hz| {
            response(
                &el,
                v.airway.oral.last_node(),
                v.airway.pharynx.first_node(),
                v.airway.glottis_port,
                hz,
            )
            .ln()
        })
    }
    let truth = [0.25, 0.18];
    let target = probe(truth[0], truth[1]);
    let error =
        |p: [f64; 12]| -> f64 { ((0..12).map(|i| (p[i] - target[i]).powi(2)).sum::<f64>()).sqrt() };
    let mut p = [-0.20, -0.15];
    let before = error(probe(p[0], p[1]));
    let mut score = before;
    for _ in 0..18 {
        let current = probe(p[0], p[1]);
        let mut deriv = [[0.0f64; 12]; 2];
        for j in 0..2 {
            let mut plus = p;
            plus[j] += 0.02;
            let mut minus = p;
            minus[j] -= 0.02;
            let a = probe(plus[0], plus[1]);
            let b = probe(minus[0], minus[1]);
            for k in 0..12 {
                deriv[j][k] = (a[k] - b[k]) / 0.04;
            }
        }
        let mut h = [[0.0f64; 2]; 2];
        let mut g = [0.0f64; 2];
        for j in 0..2 {
            for k in 0..12 {
                g[j] += (current[k] - target[k]) * deriv[j][k];
            }
            for l in 0..2 {
                for k in 0..12 {
                    h[j][l] += deriv[j][k] * deriv[l][k];
                }
            }
            h[j][j] += 0.02;
        }
        let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
        assert!(det > 1e-10, "Jacobian is not identifiable");
        let step = [
            (-g[0] * h[1][1] + g[1] * h[0][1]) / det,
            (-g[1] * h[0][0] + g[0] * h[1][0]) / det,
        ];
        let mut improved = false;
        for alpha in [1.0, 0.5, 0.25] {
            let trial = [
                (p[0] + alpha * step[0].clamp(-0.15, 0.15)).clamp(-0.7, 0.7),
                (p[1] + alpha * step[1].clamp(-0.15, 0.15)).clamp(-0.5, 0.5),
            ];
            let loss = error(probe(trial[0], trial[1]));
            if loss < score - 1e-8 {
                p = trial;
                score = loss;
                improved = true;
                break;
            }
        }
        if !improved {
            break;
        }
    }
    println!(
        "two-parameter mouth inversion: true={truth:?}, recovered={p:?}, loss {before:.6}->{score:.6}"
    );
    assert!(score < before * 0.05);
    assert!((p[0] - truth[0]).abs() < 0.04);
    assert!((p[1] - truth[1]).abs() < 0.04);
}
