# Ready-to-execute: ACS RightsLink permission for the Friden 2009 Kp,uu table (Route 2)

Target dataset: Friden M. et al. 2009, "Structure-Brain Exposure Relationships in Rat and Human
Using a Novel Data Set of Unbound Drug Concentrations in Brain Interstitial and Cerebrospinal
Fluids," J Med Chem 52(20):6233-6243, DOI 10.1021/jm901036q. ~41 marketed drugs with Kp,uu,brain -
the seminal public Kp,uu set, structures (named drugs) included.

TWO legitimate paths (do both):

A. USE THE VALUES AS FACTS NOW (no permission needed). The numeric Kp,uu values for the ~41 NAMED
   marketed drugs are facts (not copyrightable); you may extract them, look up SMILES yourself from
   PubChem, train on them, and cite the paper. Do NOT re-host the ACS SI table verbatim in the
   public repo. (This is the same basis used for the 10-drug anchor set in scripts/119.)

B. PERMISSION TO REDISTRIBUTE THE TABLE (if you want to vendor it openly):
   1. Go to the article page: https://pubs.acs.org/doi/10.1021/jm901036q
   2. Click "Rights & Permissions" (or "Reprints & Permissions") -> opens Copyright Clearance
      Center RightsLink.
   3. Choose requestor type: Academic / non-profit. Intended use: "Reuse in a Journal/Magazine" or
      "Reuse in a Thesis/Dissertation" or "Reuse in a Database/Website" as applicable.
   4. Specify: "Supporting Information data table (Kp,uu,brain values), for an open-source
      non-commercial research model; redistribution under CC-BY with attribution."
   5. Portion: "Table / data set." Format: "Electronic." Submit.
   6. Academic data-table reuse is frequently granted at no charge; if RightsLink quotes a fee or
      blocks redistribution, fall back to path A (facts) or email the corresponding author.
   Reference: https://pubs.acs.org/page/copyright/permissions.html

Deliverable once obtained: add the rows to data/raw/kpuu_train.csv (smiles, kpuu, source, license)
and retrain (see kpuu_data_acquisition_guide.md). Record the granted license string in the CSV.
