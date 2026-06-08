# Persistence-target DTI module - held-out demonstration

End-to-end test of `engine/persistence_dti.py:substrate_hypothesis` on compounds NOT used to calibrate the panel. Reproduced by `scripts/105_persistence_dti_demo.py`.

## Scorecard

- **Held-out ablative (BH3-mimetic) recovered: 0/3** (none) - generalization of the senolytic channel.
- **Flavonoid-senolytic false-durable: 0** (expected 0 - flavonoid senolytics are not BH3-mimetics; the DTI channel should be silent and the L3 Tanimoto detector handles them).
- **Non-persistence substrate leaks: 0** (expected 0 - off-substrate drugs should ABSTAIN).

## Per-compound

| compound | expected | substrate hypothesis | durable? | engaged (calibrated) | capability flags | abstained (un-calibrated) |
|---|---|---|---|---|---|---|
| obatoclax | ablative_heldout | ABSTAIN | no | - | - | - |
| gossypol | ablative_heldout | ABSTAIN | no | - | - | - |
| sabutoclax | ablative_heldout | ABSTAIN | no | - | - | HDAC6 |
| fisetin | senolytic_flavonoid | ABSTAIN | no | - | - | - |
| quercetin | senolytic_flavonoid | ABSTAIN | no | - | - | - |
| chidamide | capability_heldout | ABSTAIN | no | - | - | - |
| givinostat | capability_heldout | ABSTAIN | no | - | - | - |
| tazemetostat | capability_heldout | ABSTAIN | no | - | - | - |
| galantamine | non_persistence | ABSTAIN | no | - | - | - |
| rivastigmine | non_persistence | ABSTAIN | no | - | - | - |
| atomoxetine | non_persistence | ABSTAIN | no | - | - | - |
| citalopram | non_persistence | ABSTAIN | no | - | - | - |
| lamotrigine | non_persistence | ABSTAIN | no | - | - | - |

## Reading

Only an ABLATIVE (senolytic) engagement on a calibration-passing target promotes to durable. capability/window engagements are reported as hypotheses (flags) and never auto-promoted. Engagements on un-calibrated targets are listed but IGNORED (the `abstained` column) - the engine refuses to trust a channel MAMMAL cannot route. This is the structure-computable persistence prior the design doc deferred, now gated on measured per-target calibration.
