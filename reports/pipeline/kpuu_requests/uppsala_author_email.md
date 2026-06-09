# Ready-to-send: data-reuse request to the Uppsala Kp,uu group (Route 2)

The Uppsala group (Hammarlund-Udenaes / Irena Loryan / Markus Friden) authored most of the public
Kp,uu,brain literature and curate large internal Kp,uu tables. They are the highest-yield academic
contact for a real, SMILES-linked Kp,uu dataset under a shareable license.

WHO to email (verify current addresses on the group's Uppsala University Pharmacy / Translational
PKPD page before sending):
- Prof. Margareta Hammarlund-Udenaes (group founder)
- Assoc. Prof. Irena Loryan (corresponding author on the 2022 Kp,uu review and several datasets)
- Dr. Markus Friden (first author, Friden 2009 J Med Chem 52:6233 - the seminal Kp,uu set)

------------------------------------------------------------------------
Subject: Reuse of a Kp,uu,brain dataset for an open-source CNS free-exposure model (non-commercial)

Dear Prof. Hammarlund-Udenaes and Dr. Loryan,

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
