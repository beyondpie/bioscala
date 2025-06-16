package Peak

import GRange.GenomicRange
import SZUtils.{readTable, writeStrings2File, writeMap2File}
import Bed.GFF

/**
 * Peak from MACS2.
 * @param r
 *   GenomicRange Bed coordinate: zero-based and half-open
 * @param name
 * @param int10Neglog10qval
 * @param strand
 * @param fcAtSmt
 * @param neglog10pval
 * @param neglog10qval
 * @param relaSmt
 * @param spm
 */
case class Peak(
  r: GenomicRange, name: String, int10Neglog10qval: Int, strand: String,
  fcAtSmt: Float, neglog10pval: Float, neglog10qval: Float, relaSmt: Int,
  spm: Double
) {

  def mkString(sep: String = "\t"): String = {
    List(
        r.mkString(sep),
        name,
        int10Neglog10qval,
        strand,
        fcAtSmt,
        neglog10pval,
        neglog10qval,
        relaSmt,
        spm
    ).mkString(sep)
  }
  def getSummit: GenomicRange = {
    new GenomicRange(
        chrom = r.chrom,
        startFrom = r.startFrom + relaSmt,
        endTo = r.startFrom + relaSmt + 1
    )
  }

  /**
   * Transform to GFF for usage of superEnhancer.
   * @param fea
   * @param frame
   * @return
   *   GFF object
   */
  def toGFF(fea: String = "CRE", frame: String = "."): GFF = {
    new GFF(
        seqname = r.chrom,
        source = name,
        feature = fea,
        startFromOne = r.startFrom + 1,
        endClosed = r.endTo,
        score = spm,
        strand = strand,
        frame = frame,
        group = name
    )
  }
}

type Peaks = List[Peak]

object Peak {
  def updateSPM(r: Peaks): Peaks = {
    if (r.length < 1) {
      Nil
    } else {
      val totalScore = r.map(_.int10Neglog10qval).sum.toDouble
      r.map(x =>
        x.copy(spm =
          (1_000_000.toDouble * x.int10Neglog10qval.toDouble) / totalScore))
    }
  }
  val heads: List[String] = GenomicRange.heads ::: List(
      "name",
      "int10Neglog10qval",
      "strand",
      "foldChangeAtSummit",
      "neglog10pval",
      "neglog10qval",
      "relaSummit",
      "ScorePerMillion"
  )
  def readBed(fnm: String, maxNeglog10qval: Int = 1000,
    head: Boolean = false): Peaks = {
    val r = readTable(fnm, sep = "\t", head = head)
      .map(x =>
        new Peak(
            r = new GenomicRange(x(0), x(1).toInt, x(2).toInt),
            name = x(3),
            int10Neglog10qval = x(4).toInt.min(maxNeglog10qval),
            strand = x(5),
            fcAtSmt = x(6).toFloat,
            neglog10pval = x(7).toFloat,
            neglog10qval = x(8).toFloat,
            relaSmt = if (x.length > 9) then x(9).toInt else -1,
            spm = 0.0
        ))
    updateSPM(r)
  }

  def sortByScoreHelper(r: Peak) = {
    (r.r.chrom, r.r.startFrom, r.r.endTo, -r.neglog10qval)
  }
}

object Peaks {
  def sort(r: Peaks): Peaks = {
    r.sortBy(x => (x.r.chrom, x.r.startFrom, x.r.endTo))
  }

  def write(
    p: Peaks, fnm: String, sep: String = "\t", overwrite: Boolean = true
  ): Unit = {
    if p.length < 1 then ()
    else
      val content =
        Peak.heads.mkString(sep) :: p.map(x => x.mkString(sep))
      writeStrings2File(content, fnm, overwrite)
  }

  def writeMap(
    r: Map[String, Peaks], out: os.Path, sep: String = "\t",
    overwrite: Boolean = true, chromOrd: List[String]
  ): Unit = {
    if (r.size >= 1) {
      val content = r.map((k, v) => (k, v.map(p => p.mkString(sep))))
      writeMap2File(
          content,
          out,
          overwrite,
          ordedKeys = chromOrd,
          head = Peak.heads.mkString(sep)
      )
    }
  }
}


