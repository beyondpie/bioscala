package bioscala.LightCoord.HomerMouseGenomeAnnot

import os._
import SZUtils.ifelse
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.GenomeCoord.genomeCoordIgnoreStrandOrd
import Genome.MouseGenome

type Annot = (name: String, g: GenomeCoord, group: String)

/** Returns mouse genome annotations from Homer.
  *
  * NOTE: strand "0"="+", "1"="-" following Homer. But
  * not all the element has the concrete strand info.
  * Homer just labels them as 0. For example, CpG
  * islands and Centromeres. In intergenic.ann.txt,
  * they use "+" for strand columns.
  *
  * @param fnm
  *   Homer annotation file path (os.Path)
  * @param sep
  *   default "\t"
  */
def readGenomeElement6Cols(fnm: os.Path, sep: String = "\t"): Vector[Annot] = {
    os.read.lines
    .stream(fnm)
    .map(x => x.strip().split(sep))
    .filter(x => x.nonEmpty)
    .filter(x => x.length >= 6)
    .filter(x => MouseGenome.ordChrs.contains(x(1)))
    .map(x => {
      (
          name = x.head,
          g = (
              chrom = x(1),
              coord = (x(2).toInt, x(3).toInt),
              strand = ifelse(x(4) == "0" || x(4) == "+", "+",
                  "-")
          ),
          group = x(5)
      )
    })
    .toVector
    .sortBy(_.g)(using genomeCoordIgnoreStrandOrd)
}
