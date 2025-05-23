package GRange

trait Merge[T <: GenomicRange] {
  def mergeSorted(x: List[T]): List[T]
  def merge(x: List[T], ordChrom: List[String]): List[T] = {
    val chroms = x.map(_.chrom).distinct
    val ordedChroms = ordChrom.filter(chr => chroms.contains(chr))
    val chrom2grs = x.groupBy(_.chrom).map((chrom, grs) => {
      val ordedgrs = grs.sortBy(x => (x.startFrom, x.endTo))
      val mgrs = mergeSorted(ordedgrs)
      (chrom, mgrs)
    })
    ordedChroms.flatMap(chr => {
      chrom2grs(chr)
    })
  }
}

object BedToolGenomicRangeMerge extends Merge[GenomicRange] {
  private def mergeSortedOnce(r: List[GenomicRange]): List[GenomicRange] = {
    if (r.isEmpty) {
      List()
    } else {
      val isOvped = (r.head :: r)
        .sliding(size = 2, step = 1)
        .map(x => x.head.isOverlap4Merge(x(1)))
        .toList
      val cuts = isOvped.zipWithIndex
        .filter(!_._1)
        .map(_._2)
      ((0 :: cuts) ::: List(r.length))
        .sliding(size = 2, step = 1)
        .map(x => {
          val tmp = r.slice(x.head, x(1))
          GenomicRange(chrom = tmp.head.chrom,
            startFrom = tmp.head.startFrom, endTo = tmp.last.endTo)
        })
        .toList
    }
  }
  override def mergeSorted(r: List[GenomicRange]): List[GenomicRange] = {
    var r0 = mergeSortedOnce(r)
    var l0 = r.length
    var r1 = mergeSortedOnce(r0)
    var l1 = r1.length
    while (l0 != l1) {
      r0 = r1
      r1 = mergeSortedOnce(r1)
      l0 = r0.length
      l1 = r1.length
    }
    r1
  }
}



