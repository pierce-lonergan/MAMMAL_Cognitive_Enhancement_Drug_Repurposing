# Ready-to-send: data-reuse request to the Uppsala Kp,uu group (Route 2)

The Uppsala group (Hammarlund-Udenaes / Irena Loryan / Markus Friden) authored most of the public
Kp,uu,brain literature and curate large internal Kp,uu tables. They are the highest-yield academic
contact for a real, SMILES-linked Kp,uu dataset under a shareable license.

WHO to email (addresses verified from the official Uppsala University staff profiles, uu.se/katalog,
2026-06; use the @uu.se domain - the @farmbio.uu.se variant is legacy):
- PRIMARY: Assoc. Prof. Irena Loryan - **irena.loryan@uu.se** - active Translational PKPD group lead
  and corresponding author on the 2022 Kp,uu review and several Kp,uu datasets. Address her first.
- CC: Prof. (emeritus) Margareta Hammarlund-Udenaes - **margareta.hammarlund-udenaes@uu.se** -
  group founder; now emeritus (retired) so she may be less active, but CC her as the senior contact.
- (Dr. Markus Friden, first author of Friden 2009 J Med Chem 52:6233, moved to industry/AstraZeneca
  - no current Uppsala address; reach him via LinkedIn only if you specifically want the 2009 set,
  otherwise Loryan can point you to it.)

------------------------------------------------------------------------
Subject: Reuse of a Kp,uu,brain dataset for an open-source CNS free-exposure model (non-commercial)

Dear Dr. Loryan and Prof. Hammarlund-Udenaes,

I am developing an open-source, non-commercial machine-learning model of unbound brain exposure
(Kp,uu,brain) as part of a CNS drug-repurposing project. Our current Stage-3 model uses logBB
(B3DB, CC0) only as an honest proxy; on a 10-compound measured-Kp,uu anchor set the proxy correlates
with true Kp,uu at only Spearman ~0.5, so a real Kp,uu training set would be a major improvement.

Would you be willing to share a machine-readable table (compound SMILES or InChI + measured
Kp,uu,brain, ideally with fu,brain and Kp,brain) from your published work for model training, and to
permit redistribution under a CC-BY (or CC-BY-NC) license with full attribution? I am happy to:
- sign a Data Use Agreement / Material Transfer Agreement through my institution,
- restrict use to non-commercial research,
- cite the source dataset prominently and co-acknowledge the contribution,
- share the resulting open model and any derived analyses back with your group.

If a full table is not shareable, even the marketed-reference-compound subset (with structures) would
be very valuable. Thank you very much for considering this.

Best regards,
[Your name]
[Affiliation, department]
[email / ORCID]
------------------------------------------------------------------------

After a YES: route any DUA/MTA through your institution's research-contracts/grants office; deposit
the received table as data/raw/kpuu_train.csv (columns smiles, kpuu, source, license) and run
scripts/119 / a 111-style trainer (see kpuu_data_acquisition_guide.md).
