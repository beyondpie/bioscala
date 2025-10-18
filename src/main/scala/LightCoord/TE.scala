package bioscala.LightCoord.TE

import os._
import bioscala.LightCoord.GenomeCoord.GenomeCoord
import bioscala.LightCoord.Bed.{LightBed, LightBedElement}

type TECategory = (TEclass: String, TEfamily: String, TEsubfamily: String)

/**
  * Get TECategory infomation from Homer's TE annotation string.
  * 
  * @param s String concated by "|". For example, "L1_Mus3|LINE|L1-HOMER1",
  *   which means: class "LINE", family "L1-HOMER1", subfamily "L1_Mus3"
  */
def toTECategoryFromHomerString(s: String): TECategory = {
  val t = s.split("\\|")
  (TEclass = t(1), TEfamily = t(2), TEsubfamily = t(0))
}

def mkStringTECat(x: TECategory, sep: String): String = {
    s"${x.TEclass}$sep${x.TEfamily}$sep${x.TEsubfamily}"
}

/**
  * Load Repeats from Homer's repeat file.
  *
  * Homer's repeat file is typically located at (use mm10 as example)
  * [HomerPath]/data/genomes/mm10/mm10.repeats, which is a single file
  * that records all the TEs in a species.
  * 
  * @param f Homer's repeat file
  * @return
  */
def loadHomerRepeats(f: os.Path): LightBed = {
  os.read.lines.stream(f)
    .map(x => x.strip().split("\t"))
    .filter(x => x.nonEmpty)
    .filter(x => x.length >= 5)
    .map(x => (
      g = (
        chrom = x(1),
        coord = (startFrom = x(2).toInt, endTo = x(3).toInt),
        strand = if(x(4).toInt > 0) {"-"} else {"+"}
      ), name = x(0), score = x(4).toDouble))
    .toVector
}
