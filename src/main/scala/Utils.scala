package SZUtils

// TODO
// 1. rename SZUtils as bioscala.utils
// 2. put conversions into seperate file
// 3. put different classes of functions into different modules
// 4. add bioscala as prefix for this package

import os._
import java.io.FileNotFoundException
import java.io.File

/** Write Iterable Strings to one file
  *
  * @param content
  *   Iterable Strings
  * @param to
  *   output filename
  * @param overwrite
  *   default True, cover the original file.
  * @param head
  *   String, columns for the file if it's "" then no
  *   columns as header no new line symbol
  * Return java.nio.file.FileAlreadyExistsException if
  *  output file exists and overwrite is false.
  * TODO: fix logic when overwrite is false and file exits.
  */
def writeStrings2File(
  content: Iterable[String],
  to: String,
  overwrite: Boolean = true,
  head: String = ""
): Unit = {
  val out = os.Path(to)
  if (overwrite && os.isFile(out) && os.exists(out)) {
    println(s"${to} exist, will remove it.")
    os.remove(out)
  }
  // if (!overwrite && os.exists(out)) {
  //   Unit
  // }
  val dir = out / os.up
  if (!os.exists(dir)) {
    os.makeDir.all(dir)
  }
  if (head.nonEmpty) {
    os.write(out, s"${head}\n")
    os.write.append(out, content.mkString("\n"))
  } else {
    os.write(out, content.mkString("\n"))
  }
}

def writeMap2File[K](
  content: Map[K, List[String]],
  out: os.Path,
  overwrite: Boolean = true,
  ordedKeys: List[K],
  head: String = ""
): Unit = {
  if (overwrite && os.isFile(out) && os.exists(out)) {
    println(s"${out} exist, will remove it.")
    os.remove(out)
  }
  val dir = out / os.up
  if (!os.exists(dir)) {
    os.makeDir.all(dir)
  }
  if (head.nonEmpty) {
    os.write(out, s"${head}\n")
  }

  ordedKeys.foreach { x =>
    {
      if (content.isDefinedAt(x)) {
        // println(s"save ${x}'s values to ${out}.")
        os.write.append(out, content(x).mkString("\n"))
        os.write.append(out, "\n")
      }
    }
  }
  // println(s"write to ${out}: done.")
}

def readHead(fnm: String): String = {
  os.read.lines.stream(os.Path(fnm)).slice(0, 2).head
}

/** Read all the lines from a file At most 2147483647
  * (Int.MaxValue) lines. TODO:
  *   - use syntax that more robust to skip head
  *   - use Vector or Iterable type as return.
  *
  * @param fnm
  * @param sep
  * @param skipEmptyField
  * @param head
  * @return
  */
def readTable(fnm: String, sep: String = ",",
  skipEmptyField: Boolean = false,
  head: Boolean = false): List[List[String]] = {
  val path = os.Path(fnm)
  if (!os.exists(path)) {
    println(s"file not found: ${fnm}, return empty List.")
    List()
  } else {
    os.read.lines
      .stream(path)
      .slice(if (head) 1 else 0, Int.MaxValue)
      .map(x => {
        x.split(sep)
          .toList
          .map(x => x.strip())
          .map(x => x.replaceAll("\"", ""))
          // skip empty field
          .filter(x => skipEmptyField || x.nonEmpty)
      })
      .filter(x => x.nonEmpty)
      .toList
  }
}

def readLines(
  fnm: String,
  head: Boolean = false
): List[String] = {
  val path = os.Path(fnm)
  if (!os.exists(path)) {
    println(s"file not found: ${fnm}, return empty List.")
    List()
  } else {
    os.read.lines
      .stream(path)
      .slice(if (head) 1 else 0, Int.MaxValue)
      .toList
  }
}
def hereByGit: Option[String] = {
  var r: Option[String] = None
  var root              = os.pwd
  var stop: Boolean     = false
  while ((root.toString != "/") && !stop) {
    if (os.exists(root / ".git")) {
      r = Some(root.toString)
      stop = true
    }
    root = root / up
  }
  r
}

def getFileSize(fnm: String): Long = {
  File(fnm).length
}

def mergeListOfMap[K, V](
  a: List[Map[K, List[V]]],
  b: List[K]
): Map[K, List[V]] = {
  if (b.length < 1) {
    println("Warning: keys are empty. Map() returned.")
    Map()
  } else {
    val aa = a.filter(_.size > 0)
    if (aa.length < 1) {
      println("Warning: all Maps all empty. Map() returned.")
      Map()
    } else {
      b.map(k => {
        val v = a.map { m =>
          if m.isDefinedAt(k) then m(k)
          else List()
        }.flatten
        k -> v
      }).toMap
    }
  }
}

trait TaskElement {
  val flagfnm: String
  val skip: Boolean
  def runCore(): Unit
  def run(): Unit =
    if (os.exists(os.Path(flagfnm)) & skip) then ()
    else
      runCore()
  os.proc("touch", flagfnm).call(check = false)
}

def showProgressMessage(name: String)(f: () => Unit): Unit = {
  println(f"$name begins.")
  f()
  println(f"$name finishes.")
}

class SimpleCommandTask(
  val commands: List[String],
  val logfnm: String,
  val flagfnm: String,
  val skip: Boolean,
  val check: Boolean = true,
  val verbose: Boolean = true,
  val name: String = ""
) extends TaskElement {
  val needRun            = !os.exists(os.Path(flagfnm)) || (!skip)
  val needPrint: Boolean = verbose && (name.length() > 0)
  def runCore_(): Unit   = {
    val cs = commands.map(x => os.Shellable(List(x)))
    os.proc(cs*)
      .call(
          check = check,
          stdout = os.Path(logfnm),
          stderr = os.Path(logfnm)
      )
  }

  def runCore(): Unit = needRun match {
    case true =>
      needPrint match {
        case true  => showProgressMessage(name)(runCore_)
        case false => runCore_()
      }
    case false =>
      needPrint match {
        case true  => println(s"Skip $name job.")
        case false => ()
      }
  } // end of runCore
}   // end of class SimpleCommandTask

given str2path: Conversion[String, os.Path] with {
  def apply(x: String): os.Path = os.Path(x)
}
given path2str: Conversion[os.Path, String] with {
  def apply(x: os.Path): String = x.toString
}

given fromString2PathRedirect: Conversion[String, os.PathRedirect]
with {
  def apply(x: String): os.PathRedirect = {
    PathRedirect(os.Path(x))
  }
}

given fromString2Option: Conversion[String, Option[String]] with {
  def apply(x: String): Option[String] = {
    if (x.length > 0) {
      Some(x)
    } else {
      None
    }
  }
}


object Extensions {
  extension (a: os.Path) {
    def +(b: String): String = a.toString + b
  }
}

/** Transform color hex code to RGB. RGB seperated by
  * comma. e.g., 255,0,0.
  *
  * @param hex
  *   e.g., #FF00FF. Assume "#" included.
  * @return
  */
def colorHex2RGB(hex: String): String = {
  val n = hex.substring(1)
  hex
    .substring(1)
    .sliding(2, 2)
    .map(x => Integer.parseInt(x, 16))
    .mkString(sep = ",")
}

def toStringOption(s: String): Option[String] = s match {
  case "NA" => None
  case _    => Some(s)
}

def getStatistic(values: Seq[Double], statistic: String = "mean",
  emptyValue: Double = 0.0) = {
  statistic.toLowerCase match {
    case "mean" =>
      if (values.nonEmpty) values.sum / values.size
      else emptyValue
    case "max" => if (values.nonEmpty) values.max else emptyValue
    case "min" => if (values.nonEmpty) values.min else emptyValue
    case _     =>
      throw new IllegalArgumentException(
          s"Unsupported statistic: $statistic")
  }
}

/** A wrapper of if-else to get data of same type from
  * if-else expressions.
  *
  * @param test
  * @param yesValue
  * @param noValue
  * @return
  */
def ifelse[T](test: Boolean, yesValue: T, noValue: T): T = {
  if (test) {
    yesValue
  } else {
    noValue
  }
}
