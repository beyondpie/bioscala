package SimpleBigWig

import os._
import org.broad.igv.bbfile.BBFileReader
import scala.collection.mutable.ListBuffer
import GRange.{
  GenomicRange, mouseGenomicRangeOrd}
import SZUtils.writeStrings2File
import SZUtils.getStatistic
import bioscala.LightCoord.findOvlpOneChrSorted

case class BedGraphElement(g: GenomicRange, s: Double) {
  def mkString(sep: String = "\t"): String = {
    g.mkString(sep) + sep + s.toString
  }

  override def toString(): String = {
    g.toString + "_" + s.toString
  }
}

object BedGraphElement {
  val heads = GenomicRange.heads :+ "score"
  def mkHead(sep:String = "\t"): String = heads.mkString(sep)
}

object BedGraph {
  def sortBedGraph(x: Seq[BedGraphElement]): Seq[BedGraphElement] = {
    x.sortBy(_.g)
  }

  def toBedGraph(x: Seq[BedGraphElement], outf: String): Unit = {
    writeStrings2File(
        content = x.map(i => i.toString),
        to = outf,
        head = ""
    )
  }
  // TODO: rename it as bedgraph2bigwig
  def toBigWig(bg2bw: String, chromSizef: String, bgf: String,
    bwf: String): Unit = {
    val cs =
      List(bg2bw, bgf, chromSizef, bwf).map(x => os.Shellable(List(x)))
    os.proc(cs*).call(check = true)
  }
}

/**
 * Load bigwig data for a list of chromosomes.
 * If chr has no data, it will then be removed.
 * TODO:
 * - remove chromsomes param.
 * - rename as load bigwig
 * NOTE: during mapping the bigwigs data, it may show
 *       Null PointerException on some chrs. This might due to
 *       we do not clean up chromomes with no data at all well.
 * @param bigWigFilePath
 * @param chromosome
 * @throws Exception
 */
def loadChromosomeData(bigWigFilePath: String,
  chromosomes: List[String]): Map[String, List[BedGraphElement]] = {
  // Open the BigWig file
  val bbFileReader = new BBFileReader(bigWigFilePath)

  // Verify that it's a BigWig file
  if (!bbFileReader.isBigWigFile) {
    throw new IllegalArgumentException(
        "The file is not a valid BigWig file.")
  }
  chromosomes
    .map(chr => {
      val a = loadBigWigOneChr(bigWigFilePath, chr)
      // TODO: no need to transform to List
      (chr, a.toList)
    })
    .filter(_._2.length > 0)
    .toMap
}

/**
 * Load bigwig data for a given chromosome.
 * TODO: chromosome region should be given explicitly. 
 * @param bigWigFilePath
 * @param chr
 * @throws Exception
 */
def loadBigWigOneChr(bigWigFilePath: String,
  chr: String): Vector[BedGraphElement] = {
  val bbFileReader = new BBFileReader(bigWigFilePath)
  if (!bbFileReader.isBigWigFile) {
    throw new IllegalArgumentException(
        "The file is not a valid BigWig file.")
  }
  val a = ListBuffer[BedGraphElement]()
  try {
    val wigIterator =
      bbFileReader.getBigWigIterator(chr, 0, chr, Int.MaxValue, true)
    while (wigIterator.hasNext) {
      val item = wigIterator.next()
      val bge = new BedGraphElement(
          g = new GenomicRange(chrom = chr,
              startFrom = item.getStartBase, endTo = item.getEndBase),
          s = item.getWigValue
      )
      a.addOne(bge)
    }
  } catch {
    case np: NullPointerException =>
      println(s"Null PointerException on $chr")
  }
  a.toVector
}

/**
 * Map BigWig Signals on Regions for one chromosome.
 *
 * The results will KEEP the original region order as long as keepOrder
 * is true (by default is true).
 * 
 * NOTE: this function is still slower than R's genomic ranges overlaps.
 *       - should the bigwig library support multiple regions at one time?
 *       - in general, we should have a Better implementation of Genomic ranges
 *         overlap method, like based on B-tree used in bigwig method?
 *         should check how R's GenomicRanges implement it.
 * 
 * TODO: make param Vector as Iterable for a more general usage.
 *
 * @param bw
 * @param region
 * @param statistic
 * @param emptyValue
 * @param keepOrder
 * @throws Exception
 * @return
 */
def mapBigWigOnRegionOneChr(
  bw: Vector[BedGraphElement], region: Vector[GenomicRange],
  statistic: String = "mean", emptyValue: Double = 0.0,
  keepOrder: Boolean = true
): Vector[BedGraphElement] = {
  // check conditions
  val bwchr = bw.map(_.g.chrom).distinct
  if (bwchr.length > 1) {
    throw new IllegalArgumentException(
        "Bigwig has at least two different chromosomes.")
  }
  val rchr = region.map(_.chrom).distinct
  if (rchr.length > 1) {
    throw new IllegalArgumentException(
        "Region has at least two different chromosomes.")
  }
  if (bwchr(0) != rchr(0)) {
    throw new IllegalArgumentException(
        "Bigwig and Region have different chromsomes.")
  }

  // make sure they are sorted
  val bw2 = bw.sortBy(_.g)
  val rg2 = region.sorted
  // prepare ovlp
  val s     = bw2.map(x => (x.g.startFrom, x.g.endTo))
  val q     = rg2.map(x => (x.startFrom, x.endTo))
  val ovlps = findOvlpOneChrSorted(query = q, subject = s)
  val scores = if (ovlps.length < 1) {
    Vector.fill[Double](q.length)(emptyValue)
  } else {
    val qId2sIds = ovlps.toMap
    q.indices
      .map(i => {
        if (qId2sIds.contains(i)) {
          val sl = qId2sIds(i).startFrom
          val sr = qId2sIds(i).endTo
          getStatistic(Seq.range(sl, sr).map(i => bw2(i).s), statistic,
              emptyValue)
        } else { 0.0 }
      })
  }
  val bgs = rg2.zip(scores).map((g, s) => new BedGraphElement(g, s))
  if (!keepOrder) {
    bgs
  } else {
    // force the order to match the input region
    val rmap = bgs.map(x => (x.g, x)).toMap
    region.map(r => rmap(r))
  }
}

/**
 * Map BigWig Signals on Region for multiple chromosomes.
 *
 *   1. If some chrs in regions are not in bigwigs, they will be
 *      ignored.
 *      TODO: let chrs not in bigwigs as default emptyValues. 
 *   2. WITHIN ONE CHR, the region orders are KEPT as long
 *      as keepOrder is true (by default it's true).
 * @param bw
 * @param region
 * @param statistic
 * @param emptyValue
 * @return
 */
def mapBigWigOnRegion(bw: Map[String, Vector[BedGraphElement]],
  region: Map[String, Vector[GenomicRange]], statistic: String = "mean",
  emptyValue: Double = 0.0,
  keepOrder: Boolean = true): Map[String, Vector[BedGraphElement]] = {

  val chrs = bw.keySet ++ region.keySet
  val rg2 = region
    .filter((chr, _) => chrs.contains(chr))
  val bw2 = bw
    .filter((chr, _) => chrs.contains(chr))
  chrs.toVector
    .map(chr => {
      (chr,
          mapBigWigOnRegionOneChr(bw2(chr), rg2(chr), statistic,
              emptyValue, keepOrder))
    })
    .toMap
}

/**
  * Map BigWig on arbitrary genomic regions across single / multiple chrom(s).
  * - All the regions will NOT be ignored.
  * @param bw
  * @param region
  * @param statistic
  * @param emptyValue
  * @param keepOrder
  * @return
  */
def mapBigWigOnRegion(
  bw: Map[String, Vector[BedGraphElement]],
  region: Iterable[GenomicRange],
  statistic: String, emptyValue: Double,
  keepOrder: Boolean): Vector[BedGraphElement] = {
  val ordRegion = region.groupBy(_.chrom)
    .map((k, v) => (k, v.toVector.sorted(using mouseGenomicRangeOrd)))
  val scores =  ordRegion.map((k,v) => {
    if (!bw.contains(k)) {
      (k, v.map(x => new BedGraphElement(g = x, s = emptyValue)))
    } else {
      (k, mapBigWigOnRegionOneChr(bw(k),
        v, statistic, emptyValue, keepOrder))
    }
  })
  val r =  scores.values.flatten.toVector
  if (keepOrder) {
    val region2score = r.map(x => (x.g, x.s)).toMap
    region
      .map(x => new BedGraphElement(x, region2score(x)))
      .toVector
  } else {
    r
  }
}

/**
 * Calculates a statistic for a given genomic region from the cached
 * chromosome data. FIXME: Currently, we assume it follows bed format,
 * but bigwig can have different based, which can be got from header.
 *
 * @param chromosome
 *   Chromosome name
 * @param start
 *   Start position
 * @param end
 *   End position
 * @param mode
 *   Statistic mode: "mean", "max", "min", or "sum"
 * @return
 *   The calculated statistic
 */
def calculateBigWigRegionStatistic(
  chromosomeCache: Map[String, List[BedGraphElement]],
  chromosome: String, startFrom: Int, endTo: Int, statistic: String,
  emptyValue: Double = 0.0): Double = {
  // Check if chromosome is loaded
  if (!chromosomeCache.contains(chromosome)) {
    throw new IllegalArgumentException(
        s"Chromosome $chromosome is not loaded.")
  }

  val values: Seq[Double] = chromosomeCache(chromosome)
    .takeWhile(p => p.g.startFrom < endTo)
    .filter(p => p.g.endTo >= startFrom)
    .flatMap(item => {
      val overlapStart = item.g.startFrom.max(startFrom)
      val overlapEnd   = item.g.endTo.min(endTo)
      val v: Double    = item.s
      Seq.fill(overlapEnd - overlapStart)(v)
    })
    .toSeq
  getStatistic(values, statistic)
}
