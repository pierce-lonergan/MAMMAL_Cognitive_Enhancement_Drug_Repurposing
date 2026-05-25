"""Tests for the UniProt + PubChem + ChEMBL fetchers, using respx to mock HTTP."""

from __future__ import annotations

import httpx
import pytest
import respx

from mammal_repurposing.fetchers.chembl import top_binders, uniprot_to_chembl_target
from mammal_repurposing.fetchers.pubchem import fetch_smiles
from mammal_repurposing.fetchers.uniprot import (
    UniprotFetchError,
    fetch_sequence,
)


# --- UniProt ----------------------------------------------------------------

@respx.mock
def test_uniprot_fetch_sequence_happy_path():
    payload = {
        "primaryAccession": "P36544",
        "sequence": {"value": "MRGTPLLLVVSLFSLLQDIA", "length": 20},
        "genes": [{"geneName": {"value": "CHRNA7"}}],
        "uniProtKBCrossReferences": [
            {"database": "Ensembl",
             "properties": [{"key": "GeneId", "value": "ENSG00000175344"}]},
        ],
    }
    respx.get("https://rest.uniprot.org/uniprotkb/P36544.json").mock(
        return_value=httpx.Response(200, json=payload),
    )

    entry = fetch_sequence("P36544")
    assert entry["sequence"] == "MRGTPLLLVVSLFSLLQDIA"
    assert entry["length"] == 20
    assert entry["gene_name"] == "CHRNA7"
    assert entry["ensembl_gene_id"] == "ENSG00000175344"


@respx.mock
def test_uniprot_fetch_404_raises():
    respx.get("https://rest.uniprot.org/uniprotkb/BOGUS999.json").mock(
        return_value=httpx.Response(404),
    )
    with pytest.raises(UniprotFetchError):
        fetch_sequence("BOGUS999")


@respx.mock
def test_uniprot_missing_sequence_raises():
    payload = {"primaryAccession": "P99999", "sequence": {}}
    respx.get("https://rest.uniprot.org/uniprotkb/P99999.json").mock(
        return_value=httpx.Response(200, json=payload),
    )
    with pytest.raises(UniprotFetchError):
        fetch_sequence("P99999")


# --- PubChem ----------------------------------------------------------------

@respx.mock
def test_pubchem_canonical_smiles():
    base = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/donepezil"
        "/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    )
    respx.get(base).mock(
        return_value=httpx.Response(200, json={
            "PropertyTable": {"Properties": [{
                "CID": 3152, "CanonicalSMILES": "O=C1CC2=CC=C(OC)C=C2CC1CC1CCN(CC2=CC=CC=C2)CC1"
            }]},
        }),
    )

    hit = fetch_smiles("donepezil")
    assert hit["smiles"] is not None
    assert hit["smiles_kind"] == "canonical"
    assert hit["cid"] == 3152


@respx.mock
def test_pubchem_falls_back_to_connectivity_smiles():
    """PubChem sometimes omits CanonicalSMILES but provides ConnectivitySMILES."""
    base = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/mystery"
        "/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    )
    respx.get(base).mock(
        return_value=httpx.Response(200, json={
            "PropertyTable": {"Properties": [{
                "CID": 99, "ConnectivitySMILES": "CCO"
            }]},
        }),
    )

    hit = fetch_smiles("mystery")
    assert hit["smiles"] == "CCO"
    assert hit["smiles_kind"] == "connectivity"


@respx.mock
def test_pubchem_falls_through_alt_names():
    """Primary name 404s but the first alt resolves."""
    primary = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/bpn14770"
        "/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    )
    alt = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/zatolmilast"
        "/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    )
    respx.get(primary).mock(return_value=httpx.Response(404))
    respx.get(alt).mock(
        return_value=httpx.Response(200, json={
            "PropertyTable": {"Properties": [{
                "CID": 12345, "CanonicalSMILES": "CCC"
            }]},
        }),
    )

    hit = fetch_smiles("bpn14770", alt_names=["zatolmilast"])
    assert hit["smiles"] == "CCC"
    assert hit["name_queried"] == "zatolmilast"


@respx.mock
def test_pubchem_url_encodes_special_chars():
    """Names with commas (e.g. '7,8-dihydroxyflavone') must be URL-encoded."""
    # %2C is the encoded comma
    encoded_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/7%2C8-dihydroxyflavone"
        "/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    )
    respx.get(encoded_url).mock(
        return_value=httpx.Response(200, json={
            "PropertyTable": {"Properties": [{
                "CID": 1880, "CanonicalSMILES": "Oc1ccc(-c2cc(=O)c3ccccc3o2)c(O)c1"
            }]},
        }),
    )

    hit = fetch_smiles("7,8-dihydroxyflavone")
    assert hit["smiles"] is not None


@respx.mock
def test_pubchem_unresolvable_returns_none():
    base = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/totallymadeup"
        "/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    )
    respx.get(base).mock(return_value=httpx.Response(404))

    hit = fetch_smiles("totallymadeup")
    assert hit["smiles"] is None
    assert hit["cid"] is None


# --- ChEMBL -----------------------------------------------------------------

@respx.mock
def test_chembl_uniprot_to_target():
    respx.get("https://www.ebi.ac.uk/chembl/api/data/target.json").mock(
        return_value=httpx.Response(200, json={
            "targets": [{"target_chembl_id": "CHEMBL220", "pref_name": "Acetylcholinesterase"}],
        }),
    )
    assert uniprot_to_chembl_target("P22303") == "CHEMBL220"


@respx.mock
def test_chembl_top_binders_dedupes_and_truncates():
    activities = [
        {"molecule_chembl_id": "CHEMBL1", "molecule_pref_name": "drug1",
         "canonical_smiles": "C", "standard_type": "Ki", "standard_value": "1.0"},
        # duplicate molecule (should be deduped)
        {"molecule_chembl_id": "CHEMBL1", "molecule_pref_name": "drug1",
         "canonical_smiles": "C", "standard_type": "IC50", "standard_value": "2.0"},
        {"molecule_chembl_id": "CHEMBL2", "molecule_pref_name": "drug2",
         "canonical_smiles": "CC", "standard_type": "Ki", "standard_value": "5.0"},
        {"molecule_chembl_id": "CHEMBL3", "molecule_pref_name": "drug3",
         "canonical_smiles": "CCC", "standard_type": "Kd", "standard_value": "10.0"},
    ]
    respx.get("https://www.ebi.ac.uk/chembl/api/data/activity.json").mock(
        return_value=httpx.Response(200, json={"activities": activities}),
    )

    binders = top_binders("CHEMBL220", n=2)
    assert len(binders) == 2
    ids = [b["molecule_chembl_id"] for b in binders]
    assert ids == ["CHEMBL1", "CHEMBL2"]
    assert binders[0]["standard_value_nm"] == 1.0
