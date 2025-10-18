package Math

import SZUtils.ifelse
import smile.math.MathEx
import scala.math.{log, sqrt}

def log2(x: Double): Double = {
  log(x) / log(2)
}

def scale(x: Iterable[Double],
  center: Boolean = true): Vector[Double] = {
  val v: Double = MathEx.`var`(x.toArray)
  if (v == 0.0) {
    throw new RuntimeException("variance is zero.")
  }
  if (center) {
    val m: Double = MathEx.mean(x.toArray)
    x.map(y => (y - m) / sqrt(v)).toVector
  } else {
    x.map(y => y / sqrt(v)).toVector
  }

}

def minMaxScale(x: Iterable[Double], m: Double,
  M: Double): Vector[Double] = {
  x.map(y => (y - m) / (M - m)).toVector
}

def sumScale(x: Iterable[Double]): Vector[Double] = {
  val s = x.sum
  if (s == 0.0) {
    throw new RuntimeException("sum is zero.")
  }
  x.map(y => y / s).toVector
}

def relu(x: Double): Double = {
  ifelse(x < 0.0, 0.0, x)
}

def d2(x: Vector[Double], y: Vector[Double]): Double = {
  sqrt(x.zip(y).map((a, b) => (a - b) * (a - b)).sum)
}
