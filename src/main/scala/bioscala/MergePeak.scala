package MergePeak

import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import Peak.{Peak, Peaks}
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import Genome.Genome
import SZUtils.{readTable, writeStrings2File}
import SZUtils.mergeListOfMap
import MetaData.TSCCMeta
import SZUtils.writeMap2File

trait Filter {
  def filter(a: Peaks): Peaks
}

class BlacklistFilter(val bedfnm: String) extends Filter {
  val blist: List[GenomicRange] = readTable(bedfnm, sep = "\t", head = false)
    .map(l =>
      GenomicRange(
        chrom = l(0),
        startFrom = l(1).toInt,
        endTo = l(2).toInt
      )
    )
  val rblist = GenomicRange.reorgGenomicRanges(blist)
  def filter(a: Peaks): Peaks =
    a.filter(x => !x.r.isInSortedRanges(rblist))
}

class GenomeFilter(val fanm: String) extends Filter {
  private val useChrReg = "chr[1-9]|chr1[0-9]|chrX|chrY".r
  def UnWantedChroms(a: Peaks): Peaks =
    a.filter(x => useChrReg.matches(x.r.chrom))
  println(s"Loading genome: ${fanm}")
  val genome = new Genome(fanm)
  def NFilter(a: Peaks): Peaks =
    a.filter(x => !genome.containsN(x.r))
  def filter(a: Peaks) =
    NFilter(UnWantedChroms(a))
}


// TODO:
// - merge with GenomicRangeMerge, which is more general.
trait Merge {
  val chromOrd: List[String] = GenomicRange.chromOrd
  def mergeSorted(r: Peaks): Peaks
  def merge2(a: Map[String, Peaks]): Map[String, Peaks] = {
    a.par
      .map((k, ps) => {
        println(s"merge ${k} ... ")
        val r = (
          k,
          if (ps.nonEmpty) {
            println("sorted peaks within chrom then run merge.")
            mergeSorted(ps.sortBy(x => (x.r.startFrom, x.r.endTo)))
          } else {
            List()
          }
        )
        println(s"merge ${k} is done.")
        r
      })
      .toList
      .toMap
  }

  /** Remove the overlapped peaks.
    *
    * @param a
    *   Peaks source, will remove peaks in a having overlaps with b. Sorted or
    *   not should be OK.
    * @param b
    *   Peaks reference, sorted or not should be OK.
    * @return
    *   Peaks left from a. Keep the orders in a.
    */
  def reduce(a: Peaks, b: Peaks): Peaks = {
    val c = GenomicRange.reorgGenomicRanges(b.map(_.r))
    a.filter(x => !x.r.isInSortedRanges4Merge(c))
  }
}

object BedToolMerge2 extends Merge {
  def combineSorted(r: Peaks): Peak = {
    r.head.copy(
      new GenomicRange(
        chrom = r.head.r.chrom,
        startFrom = r.head.r.startFrom,
        // important: last element may be involved in previous ones.
        endTo = r.map(_.r.endTo).max
      ),
      int10Neglog10qval = r.map(_.int10Neglog10qval).max,
      fcAtSmt = r.map(_.fcAtSmt).max,
      neglog10pval = r.map(_.neglog10pval).max,
      neglog10qval = r.map(_.neglog10qval).max,
      spm = r.map(_.spm).max
    )
  }
  private def mergeSorted_(r: Peaks): Peaks = {
    if (r.length < 1) {
      List()
    } else {
      val isOvped = (r.head :: r)
        .sliding(size = 2, step = 1)
        .map(x => x.head.r.isOverlap4Merge(x(1).r))
        .toList
      val cuts = isOvped.zipWithIndex
        .filter(!_._1)
        .map(_._2)

      ((0 :: cuts) ::: List(r.length))
        .sliding(size = 2, step = 1)
        .map(x => {
          combineSorted(r.slice(x.head, x(1)))
        })
        .toList
    }
  }
  // twice is not enough
  // a: chr1: 1-10 b: chr1:2-3, c: chr1:4-5, d: chr1:7-11
  // in the above case, if we run mergeSorted_ once
  // we would have chr1:1-10, and chr1:4-5, chr1:7-11 since
  // we only compare adjacent two elements at one time.
  // Then at second time, we will have,
  // chr1:1-10, chr1:4-5; chr1: 7-11
  // Only run third the time here, we can have
  // chr1: 1-11 finally
  def mergeSorted(r: Peaks): Peaks = {
    var r0 = mergeSorted_(r)
    var l0 = r.length
    var r1 = mergeSorted_(r0)
    var l1 = r1.length
    while (l0 != l1) {
      r0 = r1
      r1 = mergeSorted_(r1)
      l0 = r0.length
      l1 = r1.length
    }
    r1
  }
} // end of BedToolMerge2

object BestSPMerge2 extends Merge {

  /** Reduce Peaks in the same chroms
    *
    * @param a
    * @param b
    *   sorted Peaks
    * @return
    *   Peaks in a but not in b, keep order in a.
    */
  def reduce_(a: Peaks, b: Peaks): Peaks = {
    val c = b.map(_.r)
    a.filter(x => !x.r.isInSortedRanges4Merge(c))
  }

  /** Sequentially label overlap for Peaks
    *
    * @param r
    *   Peaks cannot be empty, and restricted to the same chrom, and should be
    *   sorted based on chrom position.
    * @return
    *   List of overlap status
    */
  def isOvlp(r: Peaks): List[Boolean] = {
    var context: GenomicRange = r(0).r
    r.sliding(size = 2, step = 1)
      .map(x => {
        val t1 = x(0).r.isOverlap4Merge(x(1).r)
        val t2 = x(1).r.isOverlap4Merge(context)
        if (t1 || t2) {
          // update context range
          context = GenomicRange(
            context.chrom,
            context.startFrom,
            context.endTo.max(x(1).r.endTo)
          )
          true
        } else {
          // reset context to next one
          context = x(1).r
          false
        }
      })
      .toList
  } // end of function Ovlp

  /** One-round merge sorted Peaks in the same chrom
    *
    * @param r
    *   cannot be empty
    * @return
    *   merged peaks
    */
  def mergeSorted_(r: Peaks): Peaks = {
    val isOvped: List[Boolean] = isOvlp(r(0) :: r)
    val cuts: List[Int] = isOvped.zipWithIndex
      .filter(!_._1)
      .map(_._2)
    ((0 :: cuts) ::: List(r.length))
      .sliding(size = 2, step = 1)
      .map(x => {
        r.slice(x(0), x(1))
          .maxBy(_.spm)
      })
      .toList
  }

  def mergeSorted(r: Peaks): Peaks = {
    var r0: Peaks = mergeSorted_(r)
    // reduce keep r1 sorted
    var r1: Peaks = reduce_(r, r0)
    while (r1.length > 0) {
      val r2 = mergeSorted_(r1)
      r0 = (r0 ::: r2).sortBy(x => (x.r.startFrom, x.r.endTo))
      r1 = reduce_(r1, r0)
    }
    r0
  }
} // end of BestSPMerge2

def readReproPeak(fnm: String, head: Boolean = true): Peaks =
  readTable(fnm, sep = "\t", head = head)
    .map(x =>
      new Peak(
        r = new GenomicRange(x(0), x(1).toInt, x(2).toInt),
        name = x(3),
        int10Neglog10qval = x(4).toInt,
        strand = x(5),
        fcAtSmt = x(6).toFloat,
        neglog10pval = x(7).toFloat,
        neglog10qval = x(8).toFloat,
        relaSmt = x(9).toInt,
        spm = x(10).toDouble
      )
    )

extension (a: GenomicRange) {
  def existsInSortedRanges(c: Map[String, List[GenomicRange]]): Boolean = {
    if (c.isDefinedAt(a.chrom)) {
      val d = c(a.chrom)
      d.find(x =>
        (x.chrom == a.chrom) && (x.startFrom == a.startFrom)
          && (x.endTo == a.endTo)
      ).isDefined
    } else {
      false
    }
  }
} // end of extension for GenomicRange

object TestBestSPMerge {
  def map(a: Peaks): Map[String, Peaks] = {
    a.groupBy(_.r.chrom)
      .map((k, ps) => (k, ps.sortBy(x => (x.r.startFrom, x.r.endTo))))
  }
  val chromOrd = GenomicRange.chromOrd
  val projd = TSCCMeta.projd
  val workd = s"${projd}/04.peakcalling"
  val rscd = s"${projd}/src/test/resource"

  val sc1 = "001_CLA_EPd_CTX_Car3_Glut-H3K4me1-Male"
  val sc2 = "002_IT_EP_CLA_Glut-H3K4me1-Male"
  val sc3 = "003_L5_6_IT_TPE_ENT_Glut-H3K4me1-Male"
  val mp1 = map(readReproPeak(s"${rscd}/${sc1}.reproPeak"))
  val mp2 = map(readReproPeak(s"${rscd}/${sc2}.reproPeak"))
  val mp3 = map(readReproPeak(s"${rscd}/${sc3}.reproPeak"))

  val mp = mergeListOfMap(List(mp1, mp2, mp3), chromOrd)
  val mpp =
    mp.map((k, v) => (k, v.sortBy(_.r)))
      .map((k, v) => (k, v.map(x => x.mkString("\t"))))
      .toList
      .toMap
  SZUtils.writeMap2File(
    mpp,
    out = os.Path(s"${rscd}/allPeaksBeforeMerge.tsv"),
    ordedKeys = chromOrd,
    head = "hello"
  )
  mp.map((k, v) => (k, v.length))
  mp.map((k, v) => v.length).sum

  // Test BestSPMerge2
  val p = BestSPMerge2.merge2(mp)
  val p2 = BestSPMerge2.merge2(p)
  val p3 = BestSPMerge2.merge2(p2)
  val p4 = BestSPMerge2.merge2(p3)
  p.map((k, v) => (k, v.length))
  p.map((k, v) => v.length).sum
  p2.map((k, v) => v.length).sum
  p3.map((k, v) => v.length).sum
  p4.map((k, v) => v.length).sum
  val pp =
    p.map((k, v) => (k, v.sortBy(_.r)))
      .map((k, v) => (k, v.map(x => x.mkString("\t"))))
      .toList
      .toMap

  SZUtils.writeMap2File(
    pp,
    out = os.Path(s"${rscd}/myBestSPMerge.txt"),
    ordedKeys = chromOrd,
    head = "hello"
  )

  val rMerge =
    readTable(s"${rscd}/mergedPeaks.byR.tsv", sep = "\t", head = true)
      .map(x => GenomicRange(x(0), x(1).toInt, x(2).toInt))

  val rp = GenomicRange.reorgGenomicRanges(rMerge)

  // chr2 1331th: chr2:32170566-32172289.
  import scala.util.boundary, boundary.break
  boundary {
    for (chrom <- chromOrd) {
      val size = p(chrom).length
      for (i <- 0.until(size)) {
        val a1 = p(chrom)(i).r
        val a2 = rp(chrom)(i)
        if ((a1.startFrom != a2.startFrom) || (a1.endTo != a2.endTo)) {
          println(s"${chrom} ${i}th: ${a1}.")
          break()
        }
      }
    }
  }

  p("chr2").slice(1330, 1333)
  rp("chr2").slice(1330, 1333)

  val r = mp("chr2").sortBy(x => (x.r.startFrom, x.r.endTo))
  val content = r.map(x => x.r.toString)
  SZUtils.writeStrings2File(
    content = content,
    to = s"${rscd}/allPeaksChr2.tsv"
  )
  var r0: Peaks = BestSPMerge2.mergeSorted_(r)
  r0.map(_.r).slice(1300, 1333)
  val isOvped: List[Boolean] = BestSPMerge2.isOvlp(r(0) :: r)
  isOvped.slice(1630, 1634)
  content.slice(1630, 1634)

  val l = r.slice(1630, 1634)
  var context: GenomicRange = l(0).r
  l.sliding(size = 2, step = 1).map(x => {
    println(s"Before ${context.toString}")
    println(s"Current x(0): ${x(0).r.toString}")
    println(s"Current x(1): ${x(1).r.toString}")
    val t1 = x(0).r.isOverlap4Merge(x(1).r)
    val t2 = x(1).r.isOverlap4Merge(context)
    val b = if (t1 || t2) {
      // update context range
      context = GenomicRange(
        context.chrom,
        context.startFrom,
        x(0).r.endTo.max(x(1).r.endTo)
      )
      true
    } else {
      // reset context to next one
      context = x(1).r
      false
    }
    println(s"After ${context.toString}.")
    b
  }).toList


}
