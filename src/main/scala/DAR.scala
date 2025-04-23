package DAR
// TODO: change package name

import GenomicRange.GenomicRange


/** A simple differential genomic range.
  * @param x
  * @param log2fd
  * @param p
  * @param q
  */
case class DAR(x: GenomicRange, log2fd: Double, p: Double, q: Double)
