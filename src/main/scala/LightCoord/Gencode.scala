package bioscala.LightCoord.GenCode

import Math.log2
import SZUtils.ifelse
import java.lang.IllegalArgumentException
import bioscala.LightCoord.GenomeCoord._
import bioscala.LightCoord.GenomeCoord.GenomeCoord
import bioscala.LightCoord._

/** GFF3 line-level element uing named tuple. Ref: GTF
  * explanation
  * https://www.gencodegenes.org/pages/data_format.html
  *
  * Fields :: g: Genomic Coord including chrom,
  * startFrom, endTo, and strand gname: gene symbol
  * name gid: unique identification for this element
  * parentId:
  *   - "transcript" parent is gene
  *   - "exon" parent is transcript feaType:
  *   - feaType
  *     {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
  *   - bioType: gene_type or transcript_type such as
  *     protein_coding, pseudogene, and so on. see
  *     details here:
  *     https://www.gencodegenes.org/pages/biotypes.html
  *     geneId:
  *   - if feaType is gene, then it will be the same as
  *     gid
  *   - otherwise, will the the gene id of its parent
  *     (no matter exon or transcript)
  */
type GFF3Element = (
  g: GenomeCoord,
  gname: String,
  gid: String,
  parentId: Option[String],
  feaType: String,
  bioType: String,
  geneId: String
)

type GFF3 = Vector[GFF3Element]

def readGeneFromGenCodeGFF(gff3: String): GFF3 = {
  os.read.lines
    .stream(os.Path(gff3))
    .filter(x => !x.startsWith("#"))
    .map(x => x.strip().split("\t"))
    .filter(x => x.length >= 7)
    .map(x => {
      val meta = x.last
        .strip()
        .split(";")
        .map(x => {
          val kv = x.split("=")
          (kv(0), kv(1))
        })
        .toMap
      (
          g = (chrom = x(0),
              coord = (startFrom = x(3).toInt,
                  endTo = x(4).toInt), strand = x(6)),
          gname = meta("gene_name"),
          gid = meta("ID"),
          parentId = meta.get("Parent"),
          feaType = x(2),
          bioType = meta("gene_type"),
          geneId = meta("gene_id")
      )
    })
    .toVector
}

def getGFF3Genes(x: GFF3): GFF3 = {
  x.filter(x => x.feaType == "gene")
}

/** Get regions of genes giving upstreamd and
  * downStream size. Designed for promoter, genebody
  * and so on.
  *
  * @param g
  *   vector of GFF3Element
  * @param upStream
  *   no less than zero
  * @param downStream
  *   no less than zero
  * @param feaType
  *   name for the region, like "promoter" or "gene
  *   body"
  * @param crd3to5
  *   See details at [[LightCoord.CoordOpt]]
  * @param crd5to3
  *   See details at [[LightCoord.CoordOpt]]
  * @param chr2size
  *   context parameter to record chrom sizes.
  * @return
  *   the same order of g with modified coord
  */
def aroundGene(g: GFF3, upStream: Int, downStream: Int,
  feaType: String, crd3to5: CoordOpt, crd5to3: CoordOpt)(
  using chr2size: Map[String, Int]): GFF3 = {
  if (upStream < 0 || downStream < 0) {
    throw new IllegalArgumentException(
        "upStream or downStream should be larger than zeros.")
  }
  // TODO:
  // add Coord function for this part to support low-level operations
  g.map(x => {
    val coord = x.g.strand match {
      case "-" => {
        (
            crd3to5.getDownStream(x.g.coord, downStream),
            crd3to5
              .getUpStream(x.g.coord, upStream)
              .min(chr2size(x.g.chrom))
        ).asInstanceOf[Coord]
      }
      case _ => {
        (
            crd5to3.getUpStream(x.g.coord, upStream),
            crd5to3
              .getDownStream(x.g.coord, downStream)
              .min(chr2size(x.g.chrom))
        ).asInstanceOf[Coord]
      }
    }
    val g =
      x.g.toTuple.copy(_2 = coord).asInstanceOf[GenomeCoord]
    x.toTuple
      .copy(_1 = g, _5 = feaType)
      .asInstanceOf[GFF3Element]
  })
}

/** Get Shannon entropy of a gene from it's expression
  * across multiple groups.
  *   1. firstly transform to probs: x => x / \sum(x)
  *      2. then entropy: - \sum p * \log2(p)
  * @param x
  *   gene expression across multiple groups,
  *   non-negative.
  * @return
  */
def getShannonEntropy(x: Vector[Double],
  checkValue: Boolean = false): Double = {
  if (checkValue) {
    if (x.exists(y => y < 0.0)) {
      throw IllegalArgumentException(
          "Find negative values from the gene expressions.")
    }
    if (x.sum <= 0.0) {
      throw IllegalArgumentException(
          "No gene expressions among the gene expressions.")
    }
  }
  x.map(y => y / x.sum).map(p => -p * log2(p)).sum
}
