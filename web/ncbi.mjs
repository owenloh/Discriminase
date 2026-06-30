// NCBI E-utilities, called straight from the browser (CORS is open on eutils).
// No backend. Used to search candidate genomes and download sequences.

const EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";

function withId(url, email) {
  return url + "&tool=discriminase" + (email ? `&email=${encodeURIComponent(email)}` : "");
}
async function getJSON(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`NCBI ${r.status}`);
  return r.json();
}

// List candidate assemblies for a name -- never auto-picks; the UI shows these.
export async function searchCandidates(query, { retmax = 25, refseqOnly = true, email } = {}) {
  let term = `(${query}[Organism] OR ${query}[Title]) AND (complete genome[Title])`;
  if (refseqOnly) term += " AND (srcdb_refseq[PROP])";
  const es = await getJSON(withId(
    `${EUTILS}/esearch.fcgi?db=nucleotide&retmode=json&retmax=${retmax}&term=${encodeURIComponent(term)}`, email));
  const ids = es.esearchresult?.idlist || [];
  if (!ids.length) return [];
  const su = await getJSON(withId(
    `${EUTILS}/esummary.fcgi?db=nucleotide&retmode=json&id=${ids.join(",")}`, email));
  const out = ids.map((id) => {
    const s = su.result[id] || {};
    return {
      accession: s.accessionversion || s.caption || id,
      title: s.title || "",
      length_bp: Number(s.slen || 0),
      taxid: String(s.taxid || ""),
    };
  });
  out.sort((a, b) => b.length_bp - a.length_bp);
  return out;
}

function fastaToSeq(text) {
  return text.split("\n").filter((l) => !l.startsWith(">")).join("").replace(/\s/g, "");
}

export async function fetchSequence(accession, { email } = {}) {
  const r = await fetch(withId(
    `${EUTILS}/efetch.fcgi?db=nucleotide&id=${encodeURIComponent(accession)}&rettype=fasta&retmode=text`, email));
  if (!r.ok) throw new Error(`NCBI ${r.status} for ${accession}`);
  return { seq: fastaToSeq(await r.text()), accession };
}

export async function fetchSequenceByTaxid(taxid, { email } = {}) {
  const term = `txid${taxid}[Organism:exp] AND complete genome[Title]`;
  const es = await getJSON(withId(
    `${EUTILS}/esearch.fcgi?db=nucleotide&retmode=json&retmax=1&term=${encodeURIComponent(term)}`, email));
  const ids = es.esearchresult?.idlist || [];
  if (!ids.length) throw new Error(`no complete genome for taxid ${taxid}`);
  return fetchSequence(ids[0], { email });
}
