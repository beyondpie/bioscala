package Blacklist

import GenomicRange.GenomicRange

case class BlackList(rawBlackList: List[GenomicRange]) {
  val sortedBlackList =
    rawBlackList.sortBy(x => (x.chrom, x.startFrom, x.endTo))
}
