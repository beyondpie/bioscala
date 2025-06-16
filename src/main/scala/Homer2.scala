package Homer2

import os._
import SZUtils.TaskElement
import SZUtils.{str2path, path2str, fromString2PathRedirect}
import java.io.FileNotFoundException

/** Call motifs using HOMER2.
  *
  * homer2: http://homer.ucsd.edu/homer/homer2.html
  * http://homer.ucsd.edu/homer/introduction/install.html
  *
  * Ref:
  * http://homer.ucsd.edu/homer/ngs/peakMotifs.html.
  * NOTE: The input for bed files should be BED files
  * should have at minimum 6 columns (separated by
  * TABs, additional columns will be ignored)
  *   - Column1: chromosome
  *   - Column2: starting position
  *   - Column3: ending position
  *   - Column4: Unique Peak ID
  *   - Column5: not used
  *   - Column6: Strand (+/- or 0/1, where 0="+",
  *     1="-") In theory, HOMER will accept BED files
  *     with only 4 columns (+/- in the 4th column),
  *     and files without unique IDs, but this is NOT
  *     recommended.
  * For ChIP-Seq data, enrichment p-values typically <<
  * 1e-50.
  *
  * @param inputfnm
  *   file for motif calling (bed or other formats)
  * @param seqsize
  *   size for sequencies, -1 means "given". If
  *   analyzing ChIP-Seq peaks from a transcription
  *   factor, Chuck would recommend 50 bp for
  *   establishing the primary motif bound by a given
  *   transcription factor and 200 bp for finding both
  *   primary and "co-enriched" motifs for a
  *   transcription factor. When looking at histone
  *   marked regions, 500-1000 bp is probably a good
  *   idea (i.e. H3K4me or H3/H4 acetylated regions).
  *   In theory, HOMER can work with very large regions
  *   (i.e. 10kb), but with the larger the regions
  *   comes more sequence and longer execution time.
  *   TODO: use boolean instead of string for denovo
  * @param denovo
  *   "-nomotif" (false, don't search for de novo motif
  *   enrichment) or "" (true)
  * @param bgfnm
  *   file for background (bed or other formats). These
  *   will still be normalized for CpG% or GC% content
  *   just like randomly chosen sequences and
  *   autonormalized unless these options are turned
  *   off (i.e. "-nlen 0 -noweight"). This can be very
  *   useful since HOMER is a differential motif
  *   discovery algorithm. For example, you can give
  *   HOMER a set of peaks co-bound by another factor
  *   and compare them to the rest of the peaks. HOMER
  *   will automatically check if the background peaks
  *   overlap with the target peaks using mergePeaks,
  *   and discard overlapping regions.
  * @param mask
  *   ("-mask" or ""). Use the repeat-masked sequence.
  *   Actually, this usually doesn't matter that much.
  *   Since HOMER is a differential motif discovery
  *   algorithm, common repeats are usually in both the
  *   target and background sequences. However, it is
  *   not uncommon that a transcription factor binds to
  *   a certain class of repeats, which may cause
  *   several large stretches of similar sequence to be
  *   processed, biasing the results. Usually it's
  *   safer to go with the masked version.
  * @param hyperGenomic
  *   ("-h" or ""). By default, findMotifsGenome.pl
  *   uses the binomial distribution to score motifs.
  *   This works well when the number of background
  *   sequences greatly out number the target sequences -
  *   however, if you are using "-bg" option above, and
  *   the number of background sequences is smaller
  *   than target sequences, it is a good idea to use
  *   the hypergeometric distribution instead ("-h").
  *   FYI - The binomial is faster to compute, hence
  *   it's use for motif finding in large numbers of
  *   regions.
  * @param mcheck
  *   use a given motif database fnm to check agaist denovo motis.
  *   default None.
  * @param mknown
  *   use a given motif database fnm for enrichment. default None
  */
class MotifFinderByHomer2(
  val inputfnm: String,
  val flagfnm: String,
  val logfnm: String,
  val name: String,
  val outd: String,
  val skip: Boolean = true,
  val check: Boolean = true,
  val homer2: String =
    "/projects/ps-renlab2/szu/softwares/homer/bin",
  val seqsize: Int = -1,
  val bgfnm: Option[String] = None,
  val denovo: String = "-nomotif",
  val mask: String = "-mask",
  val hyperGenomic: String = "",
  val genome: String = "mm10",
  val mcheck: Option[String] = None,
  val mknown: Option[String] = None,
  val ncore: Int = 4
) extends TaskElement {
  val size: Int | String = if (seqsize < 0) {
    "given"
  } else { seqsize }

  val bgparam: List[String] = bgfnm match {
    case Some(x) => List("-bg", x)
    case None    => List("")
  }
  val mcheckparam: List[String] = mcheck match {
    case Some(x) => {
      if (!os.exists(x)) {
        throw new FileNotFoundException(x)
      }
      List("-mcheck", x)
    }
    case None => List("")
  }

  val mknownparam: List[String] = mknown match {
    case Some(x) => {
      if (!os.exists(x)) {
        throw new FileNotFoundException(x)
      }
      List("-mknown", x)
    }
    case None => List("")
  }
  val cs: List[String] =
    List((homer2 / "findMotifsGenome.pl").toString, inputfnm,
        genome, outd, "-size", size.toString, mask, denovo,
        hyperGenomic, "-p", ncore.toString,
        "-homer2") ++ bgparam ++ mcheckparam ++ mknownparam
  val css =
    cs
      .map(x => x.strip())
      .filter(x => x.length() > 0)
      .map(x => os.Shellable(Seq(x)))

  def runCore(): Unit = {
    println(s"Run Homer2 for $name .")
    if (!os.exists(outd)) {
      os.makeDir(outd)
    }
    os.proc(css *)
      .call(
        check = check,
        stdout = logfnm,
        stderr = logfnm
      )
  } // end of runCore
}   // end of class CallMotifByHomer2
