# CHIMERA vowel-source baseline for Aural (2026-11-09)

## Purpose

The user's imitation of Japanese「あー」sounds synthetic and rough. This test characterizes CHIMERA *before* adding another corrective effect. It does not claim improved likeness to a human reference: the user's microphone PCM was not available in the source checkout. The existing glottal model and calibrated vowel tables are **unchanged** in this experiment.

Test: `crates/chimera/tests/vowel_source_diagnostics.rs`. CHIMERA's physical tract at 48 kHz, 2× oversampling, adult body 25% female, f0 150Hz, breathiness 0.08, pressure 8, roughness 0, vibrato/drift 0, 1s rendered. FFT on steady 8192 samples, Hann window; H1−H2 from 150Hz/300Hz harmonic power. Vocal-tract F1–F3 from its **own physical transfer function**, not unstable LPC estimates on PCM.

| IPA vowel | Brightness | Tension | F1 / F2 / F3 (Hz) | H1−H2 (dB) | 2–4kHz (dB rel 0–8kHz) | 4–8kHz (dB rel 0–8kHz) |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| ɑ | .40 | .55 | 757 / 1389 / 2617 | 0.00 | −28.93 | −29.53 |
| a | .40 | .55 | 749 / 1520 / 2640 | −0.03 | −26.75 | −29.60 |
| ɑ | .00 | .55 | 757 / 1389 / 2617 | −0.44 | −34.60 | −34.06 |
| ɑ | .25 | .55 | 757 / 1389 / 2617 | −0.13 | −31.16 | −32.32 |
| ɑ | .40 | .35 | 757 / 1389 / 2617 | +2.12 | −30.99 | −31.17 |
| ɑ | .40 | .75 | 757 / 1389 / 2617 | −1.71 | −24.19 | −24.50 |

**Findings:** [a] and [ɑ] change actual vocal-tract resonances (F2 difference ~131Hz). Lowering tension increases H1−H2 and reduces higher-band energy, while reduced brightness attenuates higher-band energy. All are **controlled responses**, not evidence one setting matches the user's voice. The dependent Aural application may now expose these controls and evaluate the tone subjectively.

Japanese /a/ is generally central rather than exactly the front [a] or back [ɑ]; CHIMERA presently has both endpoints. Neither is validated as a person's Japanese vowel. The test asserts relative changes, not a universal human target.

Scientific context: [Shue et al. (2010)](https://doi.org/10.21437/Interspeech.2010-5) reports relationships between glottal open quotient, H1−H2, voice quality, and gaps. The interpretation of H1−H2 is **not independent of vocal-tract resonances**, so it should be compared at similar vowel and pitch.

## Reproduce

~~~sh
cd ~/Desktop/CHIMERA/chimera-vocal-synth
cargo test --release -p chimera --test vowel_source_diagnostics -- --nocapture
~~~

This adds one deterministic regression test to the experimental fork, without modifying the upstream source model or generated calibration tables. The wider upstream singer/speech judges were not rerun because there is no new engine sound change and the CSD/ASR datasets are not present. A previous entire-workspace test timed out on a long golden-voice run; do not treat it as a verified complete pass.
