package GRange

// case class seems to be resource cost when we have thousands of
// such elements. Now we are using named tuple (stable feature since 3.7)
// to gradually replace GenomicRanges.
// See details at LightCoord.scala

// But for small number of genomic ranges, case class is more powerful and easier
// than named tuple for usage.

// So we keep both of them at this moment.

import scala.util.matching.Regex
import scala.collection.Searching.Found
import scala.collection.Searching.InsertionPoint
import scala.math.floor
import Genome.MouseGenome
import scala.collection.mutable.ListBuffer

/**
  *  Following bed format, and no strand info implemented in this version.
  *  - index start from 0
  *  - left closed and right open 
  *
  * @param chrom
  * @param startFrom
  * @param endTo
  */
case class GenomicRange(chrom: String, startFrom: Int, endTo: Int) {
  override def toString: String =
    s"${chrom}:${startFrom}-${endTo}"

  def isOverlapFromLeft(fromLeft: GenomicRange): Boolean = {
    if (chrom == fromLeft.chrom) & (startFrom < fromLeft.endTo) then true
    else false
  }

  def isOverlapFromRight(fromRight: GenomicRange): Boolean = {
    if (chrom == fromRight.chrom) & (endTo > fromRight.startFrom) then true
    else false
  }

  def isOverlap(b: GenomicRange): Boolean = {
    if (chrom != b.chrom) {
      false
    } else if (endTo <= b.startFrom) {
      false
    } else if (startFrom >= b.endTo) {
      false
    } else {
      true
    }
  }

  def isOverlap4Merge(b: GenomicRange): Boolean = {
    if (chrom != b.chrom) {
      false
    } else if (endTo < b.startFrom) {
      false
    } else if (startFrom > b.endTo) {
      false
    } else {
      true
    }
  }

  def within(b: GenomicRange): Boolean = {
    if (!isOverlap(b)) {
      false
    } else if (startFrom >= b.startFrom && endTo <= b.endTo) {
      true
    } else {
      false
    }
  }

  def contains(b: GenomicRange): Boolean = {
    if (!isOverlap(b)) {
      false
    } else if (startFrom <= b.startFrom && endTo >= b.endTo) {
      true
    } else {
      false
    }
  }

  /**
    * Get the intersected region.
    * If no overlap, None will be returned.
    *
    * @param b
    * @return
    */
  def overlap(b: GenomicRange): Option[GenomicRange] = {
    if (!isOverlap(b)) {
      None
    } else {
      Option(GenomicRange(b.chrom, startFrom.max(b.startFrom), endTo.min(b.endTo)))
    }
  }

  /**
    * Strictly equals: all elements have to be equal.
    *
    * @param that
    * @return
    */
  override def equals(that: Any): Boolean = that match {
    case a: GenomicRange => {
      if((chrom == a.chrom) && (startFrom == a.startFrom) && (endTo == a.endTo)) {
        true
      } else {
        false
      }
    }
    case _               => false
  }

  /**
    * Compare chromosome firslty (current implementation, use String Order)
    * Then compare the coordinate (from 5' -> 3', the former the smaller).
    *
    * @param b
    * @return
    */
  def compare(b: GenomicRange): Int = {
    val t1 = chrom.compare(b.chrom)
    if t1 != 0 then t1
    else if endTo < b.startFrom then -1
    else if startFrom > b.endTo then 1
    else 0
  }

  def isInSortedRanges(c: Iterable[GenomicRange]): Boolean = {
    c.exists(b => this.isOverlap(b))
  }

  def isInSortedRanges(c: Map[String, Iterable[GenomicRange]]): Boolean = {
    if c.contains(this.chrom) then isInSortedRanges(c(this.chrom))
    else false
  }

  def isInSortedRanges4Merge(c: List[GenomicRange]): Boolean = {
    c.exists(b => this.isOverlap4Merge(b))
  }

  def isInSortedRanges4Merge(c: Map[String, List[GenomicRange]]): Boolean = {
    if c.contains(this.chrom) then isInSortedRanges4Merge(c(this.chrom))
    else false
  }

  def mkString(sep: String = "\t"): String = {
    (chrom, startFrom, endTo).toList.mkString(sep)
  }

  def isInUnWantedChrom(uwc: Regex = "chrM|chrUn|random".r): Boolean = {
    uwc.matches(chrom)
  }

  def getMidGR(): GenomicRange = {
    new GenomicRange(
      chrom = chrom,
      startFrom = getMid,
      endTo = getMid + 1
    )
  }

  def getMid: Int = {
    (startFrom + endTo) / 2
  }

  // TODO: use chrom size to check the region
  def biExtMid(extsize: Int = 2500): GenomicRange = {
    new GenomicRange(
      chrom = chrom,
      startFrom = getMid - extsize,
      endTo = getMid + extsize
    )
  }
  // def isProximal: Boolean = ???
  // def isDistal: Boolean = !isProximal

  /**
   * Get Index(0-based) for a Genomic Rnage
   *
   * @param r
   * @param leftShift  
   * @param binSize
   * @param ignoreRightIfNotFull 
   * @return
   */
  def getGenomicBinIndex(binSize: Int = 100, leftShift: Int = 0,
    ignoreRightIfNotFull: Boolean = true): Array[Int] = {
    val leftBinId: Int  = (startFrom - leftShift) / binSize
    val tmpE: Int = (endTo - leftShift) / binSize
    val rightBinId: Int = if (
      ((endTo-startFrom) % binSize > 0) & !ignoreRightIfNotFull) {
      tmpE + 1
    } else {
      tmpE
    }
    leftBinId.until(rightBinId).toArray
  }
}

object GenomicRange {
  val heads: List[String] = List("chrom", "startFrom", "endTo")
  def mkHead(sep: String = "\t"): String = heads.mkString(sep)
  val chromOrd: List[String] = 1.to(19).toList.map(i => s"chr${i}") :::
    List("chrX", "chrY")

  def fromStr(a: String, sep: String = "\t"): GenomicRange = {
    val ar = a.split(sep)
    new GenomicRange(
      chrom = ar(0),
      startFrom = ar(1).toInt,
      endTo = ar(2).toInt
    )
  }

  /**
    * Transform chr:start-end string to Genomic Rnage.
    *
    * @param a
    * @return
    */
  def fromUCSCStr(a:String): GenomicRange = {
    val a2 = a.split(":")
    val a3 = a2(1).split("-")
    new GenomicRange(a2(0), a3(0).toInt, a3(1).toInt)
  }

  /**
    * Transform a List of GenomicRanges to a chromosome
    * majored Map with sorted GenomicRanges in that chromosome.
    * TODO: use a more explicit function name
    *
    * @param a
    * @return
    */
  def reorgGenomicRanges(
    a: List[GenomicRange]
  ): Map[String, List[GenomicRange]] = {
    if (a.length < 1) then Map()
    else
      a.groupBy(_.chrom)
        .map((chrom, xs) => (chrom,
          xs.sorted(using mouseGenomicRangeOrd).toList))
  }

  def isSorted(a: Seq[GenomicRange]): Boolean = ???

} //end of GenomicRange

// FIXME: check examples of mergeSorted in Peaks.scala
// in some rare cases, may lose some peaks
@deprecated("use binary search fast but with corner bugs.")
def isInSortedRanges_(
  a: GenomicRange,
  c: Map[String, IndexedSeq[GenomicRange]]
): Boolean = {
  if c.contains(a.chrom) then
    val pool = c(a.chrom)
    val t = pool.search(a)(using mouseGenomicRangeOrd)
    t match {
      case a: Found => true
      case b: InsertionPoint => {
        val t = b.insertionPoint
        if (t == 0) {
          if a.isOverlapFromRight(pool.head) then true
          else false
        } else if (t == pool.length) {
          if a.isOverlapFromLeft(pool.last) then true
          else false
        } else {
          if a.isOverlapFromLeft(
            pool(t - 1)) | a.isOverlapFromRight(pool(t))
          then true
          else false
        }
      }
    }
  else false
}

given mouseChrOrd: Ordering[GenomicRange] with {
  private val chr2int = MouseGenome.ordChrs.zipWithIndex.toMap
  def compare(x: GenomicRange, y: GenomicRange):Int = {
    chr2int(x.chrom).compare(chr2int(y.chrom))
  }
}

given genomeCoordStartOrd: Ordering[GenomicRange] with {
  def compare(x: GenomicRange, y: GenomicRange): Int = {
    x.startFrom.compare(y.startFrom)
  }
}

given mouseGenomicRangeOrd: Ordering[GenomicRange] =
  mouseChrOrd.orElse(genomeCoordStartOrd)

def isIn(
  a: GenomicRange,
  b: List[GenomicRange]
): Boolean = {
  val c =
    b.groupBy(_.chrom)
      .map((chrom, xs) => (chrom, xs.sorted(using mouseGenomicRangeOrd).toList))
  a.isInSortedRanges(c)
}

/**
  * Filter Region A that overlap with Region B.
  *
  * @param a GenomicRanges without needing to be sorted
  * @param b GenomicRanges without needing to be sorted
  * We will organize b as Map of chrom and sorted genomic ranges.
  * @return Any element in a that overlap with b
  */
def filterByRanges(
  a: Seq[GenomicRange],
  b: Seq[GenomicRange]
): Seq[GenomicRange] = {
  val c =
    b.groupBy(_.chrom)
      .map((chrom, xs) => (chrom, xs.sorted(using mouseGenomicRangeOrd)))
  a.filter(x => x.isInSortedRanges(c))
}

