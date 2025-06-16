package Bam

import os._
import scala.collection.mutable.LinkedHashMap

import java.io.File
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMRecordIterator
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory

def getSamReader(bamfnm: String): SamReader = {
  SamReaderFactory
    .makeDefault()
    .open(new File(bamfnm))
}

def getPileUpMap(
  bamfnm: String
): LinkedHashMap[(String, Int, Boolean), Int] = {
  val f = getSamReader(bamfnm)
  val it = f.iterator()
  val r = new LinkedHashMap[(String, Int, Boolean), Int]()
  while (it.hasNext()) {
    val read = it.next
    val k = (
      read.getReferenceName(),
      read.getStart(),
      read.getReadNegativeStrandFlag()
    )
    if (r.contains(k)) {
      r += ((k, r(k) + 1))
    } else {
      r += ((k, 1))
    }
  }
  f.close()
  r
}

def removeDuplicatedPileUp(
  inbam: String,
  outbam: String,
  cutoff: Int = 10,
  skip: Boolean = true
): Unit = {
  if (os.exists(os.Path(outbam)) && skip) {
    s"${outbam} exits, and skip."
  } else {
    val f = getSamReader(inbam)
    val it = f.iterator()
    val out: SAMFileWriter = new SAMFileWriterFactory()
      .makeSAMOrBAMWriter(f.getFileHeader(), true, new File(outbam))
    // filter p by cutoff
    val p = getPileUpMap(inbam)
    p.filterInPlace((_, v) => v > cutoff)
    // output reads
    var cur_pileup = ("", 0, false)
    var count = 0
    while (it.hasNext()) {
      val read = it.next()
      val k = (
        read.getReferenceName(),
        read.getStart(),
        read.getReadNegativeStrandFlag()
      )
      if (k == cur_pileup && count != 1) {
        count -= 1
      } else if (p.contains(k)) {
        cur_pileup = k
        count = p(k)
        p -= k
      } else {
        out.addAlignment(read)
      }
    }
    out.close()
    f.close()
  }
}

def getNRead(bamfnm: String): Int = {
  var n = -1
  val reader =
    SamReaderFactory
      .makeDefault()
      .open(new File(bamfnm))
  val it = reader.iterator
  while (it.hasNext) {
    it.next
    n = n + 1
  }
  reader.close()
  n
}
