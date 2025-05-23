package bioscala.LightCoord

type LightBedElement = (g: GenomeCoord, name: String, score: Double)
type LightBed        = Vector[LightBedElement]

def mkStringLightBedElement(x: LightBedElement,
                            sep: String = "\t"): String = {
  val base = mkStringGenomeCoord(x.g, sep, useStrand = false)
  s"$base$sep${x.name}$sep${x.score.toString}$sep${x.g.strand}"
}

