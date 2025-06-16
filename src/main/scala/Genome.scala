// TODO: rename the package name
package Genome

import os._
import java.io.File
import scala.collection.mutable.HashMap
import scala.collection.mutable.ListBuffer
import scala.util.chaining._
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import SZUtils.readTable
import SZUtils.writeStrings2File

def extractChrom(line: String): String =
  line.split(" ").head.replace(">", "")

/**
 * Read Genome .fa file.
 *
 * @param fanm filenme of genome fa file.
 */
class Genome(val fanm: String) {
  val chr2seq  = HashMap[String, String]()
  val lines    = os.read.lines.stream(os.Path(fanm))
  var chr      = extractChrom(lines.head)
  var seq      = ""
  var seqlines = new ListBuffer[String]()
  for (line <- lines) {
    if (line.startsWith(">")) {
      if (seqlines.length > 0) {
        seq = seqlines.toList.mkString("")
        println(s"Add chr: ${chr} with ${seq.length} seqlines.")
        chr2seq += ((chr, seq))
        chr = extractChrom(line)
      }
      seqlines.clear()
    } else {
      seqlines += line
    }
  }

  def chrom2Size = {
    chr2seq.map((chr, seqs) => (chr, seqs.length()))
  }

  def containsN(a: GenomicRange): Boolean = {
    if chr2seq.isDefinedAt(a.chrom) then chr2seq(a.chrom)
      .slice(a.startFrom, a.endTo)
      .contains("N")
    else false
  }
}

object GenomeOpts {
  def chrom2int(x: String): Int = {
    val xchr = x.toLowerCase().replace("chr", "")
    if (xchr == "x") {
      100
    } else if (xchr == "y") {
      101
    } else {
      xchr.toInt
    }
  }
  given chromOrd: Ordering[String] with {
    def compare(x: String, y: String): Int = {
      chrom2int(x).compare(chrom2int(y))
    }
  }

}

trait GenomeImpl {
  val version: String
  val size: Long
  val ordChrs: IndexedSeq[String]
  val chr2size: Map[String, Int]
}

object MouseGenome extends GenomeImpl {
  val version        = "mm10_GRCm38_vM25"
  val mm10Size: Long = 2652783500L
  val size           = mm10Size
  val ordChrs =
    1.to(19).map(x => "chr" + x.toString) ++ List("chrX", "chrY")
  val chr2size: Map[String, Int] = Map(
      "chr1"  -> 195471971,
      "chr2"  -> 182113224,
      "chr3"  -> 160039680,
      "chr4"  -> 156508116,
      "chr5"  -> 151834684,
      "chr6"  -> 149736546,
      "chr7"  -> 145441459,
      "chr8"  -> 129401213,
      "chr9"  -> 124595110,
      "chr10" -> 130694993,
      "chr11" -> 122082543,
      "chr12" -> 120129022,
      "chr13" -> 120421639,
      "chr14" -> 124902244,
      "chr15" -> 104043685,
      "chr16" -> 98207768,
      "chr17" -> 94987271,
      "chr18" -> 90702639,
      "chr19" -> 61431566,
      "chrX"  -> 171031299,
      "chrY"  -> 91744698
  ) // end of chr2size
}

/**
 * Read mouse genome annotations from HOMER2 (Homer).
 *
 * Homer typically handle the genoem annotations under
 * [homer installed path]/data/genomes/[genome, eg mm10]/annotations
 * - basic directory: centromeres, CpG island, 3'UTR and so on.
 * - repeats:tRNA, transposase elements and so on.
 * [homer installed path]/data/genomes/mm10
 * - also has processed full annotations, repeats, tss and tts
 *
 * In this module, instead of using mm10.full.annotation, we extract
 * the corresponding annotations from the concrete files.
 */
object MouseGenomeHomerAnnot {

  /**
   * Define Annotation Atom: GenomeElement.
   *
   * @param id
   * @param range
   * @param strand "0"="+", "1"="-" following Homer. But not all the element
   *  has the concrete strand info. Homer just labels them as 0. For example,
   *  CpG islands and Centromeres.
   *  In intergenic.ann.txt, they use "+" for strand columns.
   * @param group
   */
  case class GenomeElement(id: String, range: GenomicRange,
    strand: String, group: String)

  /**
   * Group GenomeElements based on the chromosomes, then in each chrom,
   * we sorted the genome element by the start and end (increasing).
   *
   * @param x
   * @return
   */
  def groupSort(
    x: List[GenomeElement]): Map[String, List[GenomeElement]] = {
    x.groupBy(i => i.range.chrom)
      .map((k, v) => {
        (k, v.sortBy(i => i.range))
      })
  }

  def readGenomeElement6Cols(fnm: String,
    sep: String = "\t"): List[GenomeElement] = {
    readTable(fnm, sep, skipEmptyField = true, head = false)
      .map(x =>
        new GenomeElement(
            id = x(0),
            range = new GenomicRange(
                chrom = x(1),
                startFrom = x(2).toInt,
                endTo = x(3).toInt
            ),
            strand = x(4),
            group = x(5)
        ))
      .filter(x => MouseGenome.ordChrs.contains(x.range.chrom))
  }

  def readGenomeElement5Cols(fnm: String, sep: String = "\t",
    group: String): List[GenomeElement] = {
    readTable(fnm, sep, skipEmptyField = true, head = false)
      .map(x =>
        new GenomeElement(
            // add group to id in order to make it unique across multiple files.
            id = s"${group}-${x(0)}",
            range = new GenomicRange(
                chrom = x(1),
                startFrom = x(2).toInt,
                endTo = x(3).toInt
            ),
            strand = x(4),
            group = group
        ))
      .filter(x => MouseGenome.ordChrs.contains(x.range.chrom))
  }

  /**
   * Remove buffersize added to TSS or TES in Homer2,
   *  and then add new buffersize.
   *
   * @param x
   * @param bufferSize
   * @return
   */
  def reHomerBufferSize(x: GenomeElement, bufferSize: Int = 2000,
    rebufferSize: Int = 1000): GenomeElement = {
    new GenomeElement(id = x.id,
        range = new GenomicRange(
            chrom = x.range.chrom,
            startFrom =
              x.range.startFrom - 1 + bufferSize - rebufferSize,
            endTo = x.range.endTo - bufferSize + rebufferSize
        ), strand = x.strand, group = x.group)
  }

  /**
   * Map group names to unique integers to reduce memory.
   *
   * Homer uses long name like LSU-rRNA_Hsa|rRNA|rRNA-HOMER5328527 to
   * name the element (i.e., id in GenomeElement).
   *
   * @param x chromosomes and the GenomeElement
   * @return
   */
  def getAnnot2Int(
    x: Map[String, List[GenomeElement]]): Map[String, Int] = {
    val r =
      x.map((k, v) => v.map(i => i.id).distinct).toList.flatten.distinct
    val ids = 1.to(r.length)
    r.zip(ids).toMap
  }

  def getAnnot2SingleInt(v: Int = 1)(
    x: Map[String, List[GenomeElement]]): Map[String, Int] = {
    x.values.flatten
      .map(x => x.id)
      .toList
      .distinct
      .sorted
      .map(x => (x, v))
      .toMap
  }

  /**
   * Get all the genome elements from homer.
   *
   * @param homerPath directory of homer.
   * @genome default is mm10
   * @return
   */
  def getGenomeElements(homerPath: String,
    genome: String =
      "mm10"): Map[String, Map[String, List[GenomeElement]]] = {
    val genomePath = os.Path(homerPath) / "data" / "genomes" / genome
    val basicPath  = (genomePath / "annotations" / "basic").toString

    // centromere, telomere, contif, short_arm and so on.
    val GAPs =
      readGenomeElement6Cols(s"${basicPath}/gaps.ann.txt").pipe(
          groupSort)
    val CpGs =
      readGenomeElement6Cols(s"${basicPath}/cpgIsland.ann.txt").pipe(
          groupSort)
    val UTR3 =
      readGenomeElement6Cols(s"${basicPath}/utr3.ann.txt").pipe(
          groupSort)
    val UTR5 =
      readGenomeElement6Cols(s"${basicPath}/utr5.ann.txt").pipe(
          groupSort)
    // ignore UTRs, pseudogene, miRNA here.
    val Exons =
      readGenomeElement6Cols(s"${basicPath}/exons.ann.txt")
        .filter(x => x.group == "E")
        .pipe(groupSort)

    val Introns =
      readGenomeElement6Cols(s"${basicPath}/introns.ann.txt")
        .filter(x => x.group == "I")
        .pipe(groupSort)

    val TSS =
      readGenomeElement5Cols((genomePath / s"${genome}.tss").toString,
          sep = "\t", group = "TSS")
        .map(x => reHomerBufferSize(x))
        .pipe(groupSort)

    val TES =
      readGenomeElement5Cols((genomePath / s"${genome}.tts").toString,
          sep = "\t", group = "TES")
        .map(x => reHomerBufferSize(x))
        .pipe(groupSort)

    // the intergenic.ann.txt has 7 columns. but here is OK.
    // also the strand part is "+" not "0".
    val Intergenics = readGenomeElement6Cols(
        s"${basicPath}/intergenic.ann.txt").pipe(groupSort)

    Map(
        "GAP"        -> GAPs,
        "CpG"        -> CpGs,
        "3UTR"       -> UTR3,
        "5UTR"       -> UTR5,
        "Exon"       -> Exons,
        "Intron"     -> Introns,
        "TSS"        -> TSS,
        "TES"        -> TES,
        "Intergenic" -> Intergenics
    )
  }

  /**
   *  Get all the repeat elements from homer.
   *
   *   LINE.ann.txt
   *   LTR.ann.txt
   *   SINE.ann.txt
   *   Satellite.ann.txt
   *   Simple_repeat.ann.txt
   *   rRNA.ann.txt
   *   tRNA.ann.txt
   * @param homerPath
   * @param genome
   * @return
   */
  def readRepeats(homerPath: String,
    genome: String =
      "mm10"): Map[String, Map[String, List[GenomeElement]]] = {
    val genomePath = os.Path(homerPath) / "data" / "genomes" / genome
    val repeatPath = (genomePath / "annotations" / "repeats").toString
    val LINEs =
      readGenomeElement6Cols(s"${repeatPath}/LINE.ann.txt").pipe(
          groupSort)
    val LTRs =
      readGenomeElement6Cols(s"${repeatPath}/LTR.ann.txt").pipe(
          groupSort)
    val SINEs =
      readGenomeElement6Cols(s"${repeatPath}/SINE.ann.txt").pipe(
          groupSort)
    val Satellites = readGenomeElement6Cols(
        s"${repeatPath}/Satellite.ann.txt").pipe(groupSort)
    val SimpleRepeats = readGenomeElement6Cols(
        s"${repeatPath}/Simple_repeat.ann.txt").pipe(groupSort)
    val rRNAs =
      readGenomeElement6Cols(s"${repeatPath}/rRNA.ann.txt").pipe(
          groupSort)
    val tRNAs =
      readGenomeElement6Cols(s"${repeatPath}/tRNA.ann.txt").pipe(
          groupSort)
    Map(
        "LINE"         -> LINEs,
        "LTR"          -> LTRs,
        "SINE"         -> SINEs,
        "Satellite"    -> Satellites,
        "SimpleRepeat" -> SimpleRepeats,
        "rRNA"         -> rRNAs,
        "tRNA"         -> tRNAs
    )
  }

  /**
   * Map basepair to Genome Bin Id given one chromosome.
   * Start coordinate is 1.
   * @param bp
   * @param binSize
   * @return
   */
  def getBasePairId(bp: Int, binSize: Int): Int = {
    (bp / binSize) + 1
  }

  /**
   * Get the genome bin Ids for a range in one chromosome.
   * We use bed format, i.e., left closed right open for the range.
   *
   * @param startFrom
   * @param endTo
   * @param binSize
   * @return
   */
  def getGenomeBinId(startFrom: Int, endTo: Int,
    binSize: Int = 200): Vector[Int] = {
    val binIdLeft  = getBasePairId(startFrom, binSize)
    val binIdRight = getBasePairId((endTo - 1).max(startFrom), binSize)
    binIdLeft.to(binIdRight).toVector
  }

  // def getGenomeBinId(startFrom: Int, endTo: Int, binSize: Int = 200): Vector[Int] = {
  //   val binIdLeft = getBasePairId(startFrom, binSize)
  //   val binIdRight = getBasePairId(endTo, binSize)
  //   binIdLeft.to(binIdRight).toVector
  // }

  /**
   * Annotate genome bins based on the GenomeElement provided.
   *  Genome bins use bed format, i.e.,
   *  start from 0 (included), to end (excluded)
   *
   * NOTE: current implementation results a very little different with
   *  the R implementation using GenomicRnages for TSS annotation.
   *  This should be OK and will not affect the final results.
   *
   * @param chromSize for one chromosome
   * @param binSize
   * @param g sorted GenomeElements
   * @param g2i the integer index for names in GenomeElements
   *  Use the function getAnnote2Int.
   * @return
   */
  def annotateGenomeBin(chromSize: Int, binSize: Int = 200,
    g: List[GenomeElement], g2i: Map[String, Int],
    emptyId: Int = -1): IndexedSeq[Int] = {
    // toMap will unique to id (use the latest value)
    val id2gi: Map[Int, Int] = g
      .map(x =>
        getGenomeBinId(x.range.startFrom, x.range.endTo, binSize).map(
            i => (i, g2i(x.id))))
      .flatten
      .toMap

    // floor Int
    val size = (chromSize / binSize).toInt
    1.to(size)
      .map(i =>
        id2gi.contains(i) match {
          case true  => id2gi(i)
          case false => emptyId
        })
  }

  /**
   * Annote the mouse mm10 genome using homer annotaitons.
   *
   * @param args
   */
  def main(args: Array[String]) = {
    val homerPath    = "/projects/ps-renlab2/szu/softwares/homer"
    val genome       = "mm10"
    val projd        = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
    val chromSizefnm = s"${projd}/meta/mm10.chrom.sizes.lite"
    val chromSizes = readTable(chromSizefnm, "\t", head = false)
      .map(x => (x(0), x(1).toInt))
      .toMap

    val genomeElements = getGenomeElements(homerPath, genome)
    val repeatElements = readRepeats(homerPath, genome)
    val allElements    = genomeElements ++ repeatElements
    val chrs = 1.to(19).map(x => s"chr${x}") ++ List("chrX", "chrY")

    val outd      = s"${projd}/data/genome"
    val binAnnotd = s"${outd}/mm10GenomeBinAnnot"

    // for each genome element, get the g2i.
    val genomeElementKeys = genomeElements.keys.toList.sorted
    // val g2iEs = genomeElementKeys.map( x => x match {
    //   case "GAP" => ("GAP", genomeElements("GAP").pipe(getAnnot2Int))
    //   case _ => (x, genomeElements(x).pipe(getAnnot2SingleInt(1)))
    // }).toMap
    val g2iEs = genomeElementKeys
      .map(x => (x, genomeElements(x).pipe(getAnnot2SingleInt(1))))
      .toMap

    // for each repeat element, get the g2i
    val repeatKeys = repeatElements.keys.toList.sorted
    // val g2iRs = repeatKeys.map(
    //   x => (x, repeatElements(x).pipe(getAnnot2Int))
    // ).toMap
    val g2iRs = repeatKeys
      .map(x => (x, repeatElements(x).pipe(getAnnot2SingleInt(1))))
      .toMap

    // write the maps to one files (only neccessary ones)
    // Currently we use binary mapping for simplificaitons

    // annotate the genome bins for each chr, and then output to files
    val keys = genomeElementKeys ++ repeatKeys
    val g2is = g2iEs ++ g2iRs
    for (chr <- chrs) {
      println(s"Annotate chromosome ${chr}")
      val t1 = keys.map(x => {
        (x,
            annotateGenomeBin(chromSize = chromSizes(chr),
                binSize = 200, g = allElements(x)(chr), g2i = g2is(x),
                emptyId = 0))

      })
      val head = t1.map(x => x(0))
      val data = t1.map(x => x(1).toList).transpose
      writeStrings2File(
          content = data.map(x => x.mkString(",")),
          to = s"${binAnnotd}/${chr}.binary.annot.csv",
          overwrite = true,
          head = head.mkString(",")
      )
    }
  }

  def test(args: Array[String]) = {
    val homerPath    = "/projects/ps-renlab2/szu/softwares/homer"
    val genome       = "mm10"
    val projd        = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
    val chromHMMd    = s"${projd}/06.ChromHMM/out/m-bpeak-s18_pd-obs"
    val chromSizefnm = s"${projd}/meta/mm10.chrom.sizes.lite"
    val chromSizes = readTable(chromSizefnm, "\t", head = false)
      .map(x => (x(0), x(1).toInt))
      .toMap

    val genomePath = os.Path(homerPath) / "data" / "genomes" / genome
    val basicPath  = (genomePath / "annotations" / "basic").toString
    val TSS =
      readGenomeElement5Cols((genomePath / s"${genome}.tss").toString,
          sep = "\t", group = "TSS")
        .map(x => reHomerBufferSize(x))
        .pipe(groupSort)

    val chr = "chr3"
    val g   = TSS(chr)
    val g2i = TSS
      .map((k, v) => v.map(i => i.id).distinct)
      .toList
      .flatten
      .distinct
      .map(x => (x, 1))
      .toMap
    val a = annotateGenomeBin(chromSizes(chr), 200, g = g, g2i, 0)

    val sc = "001_CLA_EPd_CTX_Car3_Glut"
    val b =
      readTable(s"${chromHMMd}/${sc}_${chr}_s18-200bin-pd-obs.csv",
          sep = ",", head = true)
        .map(x => x(4))
        .map(x =>
          x match {
            case "D" => 0
            case "P" => 1
          })
    val d = a.zip(b).map((i, j) => i - j)
  }
}
