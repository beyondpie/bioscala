/**
 * Implementation of bed-related formats. Ref: bedtools
 * https://bedtools.readthedocs.io/en/latest/content/general-usage.html
 */
package Bed

import os._
import GRange.GenomicRange
import SZUtils.readTable
import SZUtils.writeStrings2File

// TODO: use g for GenomicRange in order to keep consistent with
// BedGraphElement
case class BedElement4(x: GenomicRange, name: String) {
  def mkString(sep: String = "\t"): String = {
    x.mkString(sep) + sep + name
  }
}
object BedElement4 {

  /**
   * Read bed file with four columns.
   *
   * @param fnm
   *   bed filename
   * @param head
   */
  def readBed4(fnm: String) = {
    os.read.lines
      .stream(os.Path(fnm))
      .map(x => x.strip().split("\t"))
      .map(e =>
        new BedElement4(
            new GenomicRange(e(0), e(1).toInt, e(2).toInt),
            e(3)
        ))
      .toList
  } // end of fn readBed4

  /**
   * Contruct BedElement4 from two-column file. Each line is like:
   * "chr1:3094725-3095224,Chr-O"
   * @param fnm
   * @param sep
   * @param head
   */
  def fromTwoCols(fnm: String, sep: String, head: Boolean) = {
    os.read.lines.stream(os.Path(fnm))
      .slice(if(head) 1 else 0, Int.MaxValue)
      .map(x => x.strip().split(sep))
      .map(x =>
        new BedElement4(
            GenomicRange.fromUCSCStr(x(0)),
            x(1)
        )
      ).toVector
  }
}

/**
 * Classic 10-column BEDPE format per row.
 *
 * @param x
 * @param y
 * @param name
 * @param score
 * @param strand1
 * @param strand2
 */
case class BedPE10Element(x: GenomicRange, y: GenomicRange,
  name: String = "bedpe", score: Double = 0.0, strand1: String = ".",
  strand2: String = ".") {
  def mkString(sep: String = "\t"): String = {
    List(
        x.mkString(sep),
        y.mkString(sep),
        name,
        score.toString,
        strand1,
        strand2
    ).mkString(sep)
  }
}

object BedPE10Element {
  def readBedPE(fnm: String,
    head: Boolean = false): List[BedPE10Element] = {
    readTable(fnm, sep = "\t", head = head).map(l =>
      new BedPE10Element(
          x = new GenomicRange(l(0), l(1).toInt, l(2).toInt),
          y = new GenomicRange(l(3), l(4).toInt, l(5).toInt),
          name = l(6),
          score = l(7).toDouble,
          strand1 = l(8),
          strand2 = l(9)
      ))
  } // end of fn readBedPE

  def writeBedPE(x: List[BedPE10Element], fnm: String,
    head: String = ""): Unit = {
    writeStrings2File(
        content = x.map(e => e.mkString("\t")),
        to = fnm,
        overwrite = true,
        head = head
    )
  } // end of fn writeBedPE
}

/**
 * GFF (General Feature Format) Ref:
 * https://genome.ucsc.edu/FAQ/FAQformat.html#format3 We can use it to
 * represent superEnhancer. See
 * http://younglab.wi.mit.edu/super_enhancer_code.html for details.
 * @param seqname
 * @param source
 *   unique ID
 * @param feature
 *   the type of the feature, like exons, codons
 * @param strand
 *   [+|-|.]
 * @param frame
 *   [0-2|.] reading frame of the first base
 * @param group
 *   unique ID in superEnhancer. If not equal to source, will use
 *   source.
 */
case class GFF(
  seqname: String, source: String, feature: String, startFromOne: Int,
  endClosed: Int, score: Double, strand: String, frame: String,
  group: String
) {

  /**
   * to String for writing as one line
   *
   * @param sep
   * @return
   *   a string
   */
  def toSingleStringLine(sep: String = "\t"): String = {
    List(
        seqname,
        source,
        feature,
        startFromOne.toString,
        endClosed.toString,
        score.toString,
        strand,
        frame,
        group
    ).mkString(sep)
  }
}
