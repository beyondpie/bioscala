import org.scalatest._
import org.scalatest.funsuite.AnyFunSuite
import bioscala.LightCoord.{d, Coord, getClosedCoord}

class LightCoordDistance extends AnyFunSuite {
  test("Test calculate distance") {
    assert(d((startFrom = 0, endTo = 100),
        (startFrom = 0, endTo = 10)) == 0)
    assert(d((startFrom = 100, endTo = 200),
        (startFrom = 190, endTo = 300)) == 0)
    assert(d((startFrom = 100, endTo = 200),
        (startFrom = 90, endTo = 150)) == 0)
    assert(d((startFrom = 100, endTo = 200),
        (startFrom = 40, endTo = 60)) == 40)
    assert(d((startFrom = 100, endTo = 200),
        (startFrom = 250, endTo = 300)) == 50)
  }

  test("Test getClosedCoord: no coords within range.") {
    val x = (startFrom = 500, endTo = 520)
    val y = Vector((20, 30), (40, 60), (90, 100), (550, 600)).map(
        x => x.asInstanceOf[Coord])
    assert(getClosedCoord(x, y, within = 20).isEmpty)
  }

  test("Test getClosedCoord: no overlap but within range") {
    val x = (startFrom = 500, endTo = 520)
    val y = Vector((20, 30), (40, 60), (90, 100), (550, 600),
      (100, 490)).map(
        x => x.asInstanceOf[Coord])
    val r = getClosedCoord(x, y, within = 100)
    println(r)
    assert(r.length == 1)
  }

  test("Test getClosedCoord: with overlap") {
    val x = (startFrom = 500, endTo = 520)
    val y = Vector((550, 600), (100, 600), (300, 501), (100, 200))
      .map(x => x.asInstanceOf[Coord])
    val r = getClosedCoord(x, y, within = 1000)
    println(r)
    assert(r.length == 2)
  }

}
