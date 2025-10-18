package bioscala.LightCoord.Bed

import bioscala.LightCoord.GenomeCoord._
import os._

type LightBedElement = (g: GenomeCoord, name: String,
  score: Double)
type LightBed = Vector[LightBedElement]


// TODO: use extension methods
def mkStringLightBedElement(x: LightBedElement,
  sep: String = "\t"): String = {
  val base = mkStringGenomeCoord(x.g, sep, useStrand = false)
  s"$base$sep${x.name}$sep${x.score.toString}$sep${x.g.strand}"
}


def readLightBed(fnm: String): LightBed = {
  os.read.lines
    .stream(os.Path(fnm))
    .map(x => x.strip().split("\t"))
    .filter(x => x.length == 6)
    .map(x =>
      (g = (
              chrom = x(0),
              coord = (startFrom = x(1).toInt,
                  endTo = x(2).toInt),
              strand = x.last
          ), name = x(3), score = x(4).toDouble))
    .toVector
}
