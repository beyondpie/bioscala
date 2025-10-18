package bioscala.LightCoord.BedGraph

import bioscala.LightCoord.GenomeCoord.GenomeCoord
import bioscala.LightCoord.GenomeCoord.mkStringGenomeCoord
import SZUtils.ifelse
import SZUtils.{str2path, path2str}
import java.io.FileNotFoundException
import org.broad.igv.bbfile.BBFileReader
import scala.collection.mutable.ListBuffer
import scala.util.Try
import scala.util.Success
import scala.util.Failure

type LightBedGraphElement = (g: GenomeCoord, s: Double)
type LightBedGraph        = Vector[LightBedGraphElement]

object LightBedGraphElement {
  extension (c: LightBedGraphElement) {
    def mkString(sep: String = "\t"): String = {
      val base = mkStringGenomeCoord(c.g, sep, useStrand = false)
      base + sep + c.s.toString + sep + c.g.strand
    }
  }
}

object LightBedGraph {
  def heads(useStrand: Boolean = false) = {
    val b = Vector("chrom", "startFrom", "endTo", "score")
    ifelse(useStrand, b :+ "strand", b)
  }

  /** Use ucsd bigToBigwig tool to transformer bedgraph
    * to bigwig. Default parameters will be used. Ref:
    * https://genome.ucsc.edu/goldenpath/help/bigWig.html.
    *
    * @param bg2bwBin
    * @param chromsizefnm
    * @param bgf
    * @param bwf
    */
  def toBigWig(bg2bwBin: String, chromsizefnm: String,
    bgf: String, bwf: String): Unit = {
    if (!os.exists(bgf)) {
      throw new FileNotFoundException(bgf)
    }
    if (!os.exists(chromsizefnm)) {
      throw new FileNotFoundException(chromsizefnm)
    }
    if (!os.exists(os.Path(bwf) / os.up)) {
      throw new FileNotFoundException(os.Path(bwf) / os.up)
    }
    val cs = List(bg2bwBin, bgf, chromsizefnm, bwf)
      .map(x => os.Shellable(List(x)))
    os.proc(cs*).call(check = true)
  }
}

def getBigWigReader(fnm: String): BBFileReader = {
  if (!os.exists(fnm)) {
    throw new FileNotFoundException(fnm)
  }
  val r = new BBFileReader(fnm)
  if (!r.isBigWigFile()) {
    throw new RuntimeException(s"$fnm is not a bigwig file.")
  }
  r
}

/** Returns vector of wig values and associated regions
  * in the bigwig for a genome coord.
  *
  * @param r
  *   BBFileReader See also [[getBigWigReader]]
  * @param x
  *   single genome element
  * @param emptyValue
  *   default value if wig values found
  * @param contained
  *   if false (default), any overlapped regions in the
  *   bigwig will be considered, otherwise (true) only
  *   regions within the region is considered. Use
  *   false (default) unless you have particular
  *   requests.
  */
def getWigValue(r: BBFileReader, x: GenomeCoord,
  emptyValue: Double = 0.0,
  contained: Boolean = false): LightBedGraph = {
  Try({
    val a  = ListBuffer[LightBedGraphElement]()
    val wI = r.getBigWigIterator(x.chrom, x.coord.startFrom,
      x.chrom, x.coord.endTo, contained)
    while (wI.hasNext) {
      val item = wI.next()
      a.addOne((g = (chrom = item.getChromosome(),
        coord = (startFrom = item.getStartBase(),
          endTo = item.getEndBase()),
        strand = "."),
        s = item.getWigValue()))}
    a.toVector
  }) match {
    case Success(v) => {
      if (v.nonEmpty) {
        v
      } else {
        Vector((g = x, s = emptyValue))
      }
    }
    case Failure(e) => Vector((g = x, s = emptyValue))
  }
}

/**
  * Get average value from a bigwig file for a given region.
  *
  * Here the average is the weighted sum for all the bigwig regions
  * overlapped with the given region.
  * 
  * weightByOverlapRegion: the overlapped region size / sum of all the overlapped
  * region size (NOT the size of the given region).
  * weightByGivenRegion: the overlapped region size /
  *  the size of the given region.
  * 
  * @param bwr
  * @param r
  * @param emptyValue
  * @param weightByGivenRegion
  * @param contained See [[getWigValue]]
  * @return
  */
def getbwOneRegion(
  bwr: BBFileReader, r: GenomeCoord, emptyValue: Double,
  weightByGivenRegion: Boolean, contained: Boolean = false
): Double = {
  val t = getWigValue(bwr, r, emptyValue, contained)
  if (t.length < 2) {
    t.head.s
  } else {
    val wv = t.map(
      y => {
        val width = y.g.coord.endTo.min(r.coord.endTo) -
        y.g.coord.startFrom.max(r.coord.startFrom)
        (width, y.s)
      }
    )
    val ss = wv.map(x => x._1.toDouble * x._2).sum
    if (weightByGivenRegion) {
      ss / (r.coord.endTo - r.coord.startFrom).toDouble
    } else {
      ss / wv.map(_._1).sum.toDouble
    }
  }
}
