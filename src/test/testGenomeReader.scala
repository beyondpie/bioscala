import os._
import java.io.File
import scala.collection.mutable.HashMap
import scala.collection.mutable.ListBuffer

def extractChrom(line:String):String =
  line.split(" ").head.replace(">", "")

class Genome(val fanm: String) {
  val chr2seq = HashMap[String, String]()
  val lines = os.read.lines.stream(os.Path(fanm))
  var chr = extractChrom(lines.head)
  var seq = ""
  var seqlines = new ListBuffer[String]()
  for (line <- lines) {
    if line.startsWith(">") & (seqlines.length > 0) then
      seq = seqlines.toList.mkString("")
      println(s"Add chr: ${chr} with ${seq.length} seqlines.")
      chr2seq += ((chr, seq))
      chr = extractChrom(line)
      seqlines.clear()
    else
      seqlines += line
  }
}


object TestGenomeReader{
  val proj = os.pwd
  val genomefnm = proj / "genome" / "GRCm39.genome.fa"
  val mm = new Genome(genomefnm.toString)
}
