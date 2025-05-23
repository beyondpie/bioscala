package TF

case class TranscriptionFactor (
  geneSymbol: String,
  ensembl: String,
  family: String
  )

case class TranscriptionCoFactor(
  geneSymbol: String,
  ensembl: String,
  family: String)


def readTFromAnimalTFDB(fnm: String): Vector[TranscriptionFactor] = ???
def readTCoFFromAnimalTFDB(fnm: String): Vector[TranscriptionCoFactor] = ???
