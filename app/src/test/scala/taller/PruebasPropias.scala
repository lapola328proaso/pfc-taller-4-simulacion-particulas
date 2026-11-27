package taller

import org.scalatest.funsuite.AnyFunSuite
import scala.util.Random

class PruebasPropias extends AnyFunSuite {

  val sim = new SimulacionParticulas()
  import sim._

  val eps = 1e-9
  def casi(a: Double, b: Double): Boolean = math.abs(a - b) <= eps * math.max(1.0, math.abs(b))
  def casiV(a: Vector2, b: Vector2): Boolean = casi(a.x, b.x) && casi(a.y, b.y)

  def aleatorio(n: Int, semilla: Int): Sistema = {
    val r = new Random(semilla)
    Array.fill(n)(
      Particula(
        Vector2(r.nextDouble() * 200 - 100, r.nextDouble() * 200 - 100),
        Vector2(r.nextDouble() * 20 - 10, r.nextDouble() * 20 - 10),
        r.nextDouble() * 10
      )
    )
  }
  val umbrales = List(1, 2, 3, 7, 64, 500, 1000, 5000)

  // ------------------------------ Punto 1 ------------------------------
  test("P1: dt negativo mueve hacia atrás") {
    val p = Particula(Vector2(5, 5), Vector2(1, -2), 1.0)
    assert(pasoSecuencial(Array(p), -3.0).head.pos == Vector2(2.0, 11.0))
  }
  test("P1: velocidad y masa no cambian") {
    val s = aleatorio(50, 1)
    val n = pasoSecuencial(s, 0.7)
    assert(n.map(_.vel).toList == s.map(_.vel).toList)
    assert(n.map(_.masa).toList == s.map(_.masa).toList)
  }
  test("P1: pasoSegmento devuelve solo el tramo pedido") {
    val s = Array.tabulate(10)(i => Particula(Vector2(i, 0), Vector2(1, 1), 1.0))
    val t = pasoSegmento(s, 1.0, 3, 6)
    assert(t.length == 3)
    assert(t.map(_.pos).toList == List(Vector2(4, 1), Vector2(5, 1), Vector2(6, 1)))
  }
  test("P1: dos segmentos contiguos reconstruyen el sistema (no pierde ni repite)") {
    val s = aleatorio(101, 2)
    val k = 37
    val unido = pasoSegmento(s, 0.3, 0, k) ++ pasoSegmento(s, 0.3, k, s.length)
    assert(unido.toList == pasoSecuencial(s, 0.3).toList)
  }
  test("P1: sistema vacío, secuencial y paralelo") {
    assert(pasoSecuencial(Array(), 1.0).isEmpty)
    assert(pasoParalelo(Array(), 1.0, 4).isEmpty)
  }
  test("P1: pasoParalelo == pasoSecuencial con tamaños impares y varios umbrales") {
    List(1, 2, 3, 17, 1001).foreach { n =>
      val s = aleatorio(n, n)
      umbrales.foreach(u => assert(pasoParalelo(s, 0.25, u).toList == pasoSecuencial(s, 0.25).toList))
    }
  }
  test("P1: umbral no positivo no cuelga y da lo mismo") {
    val s = aleatorio(20, 5)
    assert(pasoParalelo(s, 1.0, 0).toList == pasoSecuencial(s, 1.0).toList)
  }

  // ------------------------------ Punto 2 ------------------------------
  test("P2: energía de una partícula en reposo es 0") {
    assert(energiaSecuencial(Array(Particula(Vector2(1, 1), Vector2(0, 0), 9.0))) == 0.0)
  }
  test("P2: energía con velocidad (1,1) y masa 3 es 3") {
    assert(casi(energiaSecuencial(Array(Particula(Vector2(0, 0), Vector2(1, 1), 3.0))), 3.0))
  }
  test("P2: masaTotal suma masas y vacío da 0") {
    val s = Array(Particula(Vector2(0, 0), Vector2(0, 0), 1.5), Particula(Vector2(1, 0), Vector2(0, 0), 2.5))
    assert(masaTotal(s) == 4.0)
    assert(masaTotal(Array()) == 0.0)
  }
  test("P2: centro de masa con masas distintas pesa más a la pesada") {
    val s = Array(Particula(Vector2(0, 0), Vector2(0, 0), 1.0), Particula(Vector2(10, 0), Vector2(0, 0), 3.0))
    assert(casiV(centroDeMasaSecuencial(s), Vector2(7.5, 0.0)))
  }
  test("P2: centro de masa con masas iguales es el punto medio") {
    val s = Array(Particula(Vector2(1, 1), Vector2(0, 0), 2.0), Particula(Vector2(3, 5), Vector2(0, 0), 2.0))
    assert(casiV(centroDeMasaSecuencial(s), Vector2(2.0, 3.0)))
  }
  test("P2: masa total cero da el origen (sin NaN) y vacío también") {
    val s = Array(Particula(Vector2(4, 4), Vector2(1, 1), 0.0))
    assert(centroDeMasaSecuencial(s) == Vector2(0, 0))
    assert(centroDeMasaSecuencial(Array()) == Vector2(0, 0))
  }

  // ------------------------------ Punto 3 ------------------------------
  test("P3: energiaParalela coincide con la secuencial para todos los umbrales") {
    val s = aleatorio(1234, 3)
    umbrales.foreach(u => assert(casi(energiaParalela(s, u), energiaSecuencial(s))))
  }
  test("P3: centroDeMasaParalelo coincide con el secuencial para todos los umbrales") {
    val s = aleatorio(999, 4)
    umbrales.foreach(u => assert(casiV(centroDeMasaParalelo(s, u), centroDeMasaSecuencial(s))))
  }
  test("P3: la trampa: mitades de distinto peso no se promedian") {
    val s = Array(Particula(Vector2(0, 0), Vector2(0, 0), 1.0), Particula(Vector2(10, 0), Vector2(0, 0), 3.0))
    val r = centroDeMasaParalelo(s, 1)
    assert(casiV(r, Vector2(7.5, 0.0)))
    assert(!casiV(r, Vector2(5.0, 0.0))) // promedio de promedios, incorrecto
  }
  test("P3: sistema vacío en paralelo") {
    assert(energiaParalela(Array(), 4) == 0.0)
    assert(centroDeMasaParalelo(Array(), 4) == Vector2(0, 0))
  }
  test("P3: masa total cero en paralelo da el origen") {
    val s = Array.tabulate(9)(i => Particula(Vector2(i, i), Vector2(1, 1), 0.0))
    assert(centroDeMasaParalelo(s, 2) == Vector2(0, 0))
  }
  test("P3: tamaño impar y umbral pequeño") {
    val s = aleatorio(7, 9)
    assert(casi(energiaParalela(s, 2), energiaSecuencial(s)))
    assert(casiV(centroDeMasaParalelo(s, 2), centroDeMasaSecuencial(s)))
  }

  // ------------------------------ Punto 4 ------------------------------
  test("P4: trabajo es n") {
    assert(trabajo(0) == 0)
    assert(trabajo(5000) == 5000)
  }
  test("P4: profundidad con umbral mayor o igual a n es n") {
    assert(profundidad(10, 100) == 10)
    assert(profundidad(10, 10) == 10)
  }
  test("P4: casos pequeños y vacío") {
    assert(profundidad(0, 1) == 0)
    assert(profundidad(1, 1) == 1)
    assert(profundidad(3, 1) == 3) // 3 -> 2 -> 1: dos particiones y 1
  }
  test("P4: n no potencia de 2, sigue la mitad mayor") {
    assert(profundidad(1000, 100) == 67) // 1000,500,250,125,63: 4 particiones + 63
  }
  test("P4: tiempoEstimado con un procesador y sin partir es 2n; con 8 es n + n/8") {
    assert(tiempoEstimado(1000, 1, 1000) == 2000.0)
    assert(tiempoEstimado(1000, 8, 1000) == 1125.0)
  }
  test("P4: más procesadores nunca aumentan el tiempo estimado") {
    val ts = List(1, 2, 4, 8, 16).map(p => tiempoEstimado(1 << 20, p, 1000))
    assert(ts.zip(ts.tail).forall { case (a, b) => b < a })
  }

  // ------------------------------ Punto 5 ------------------------------
  test("P5: aceleración 5 y 0.5 (más lento)") {
    assert(aceleracion(50.0, 10.0) == 5.0)
    assert(aceleracion(30.0, 60.0) == 0.5)
  }
  test("P5: tiempo secuencial 0 da aceleración 0") {
    assert(aceleracion(0.0, 5.0) == 0.0)
  }
  test("P5: lista vacía devuelve el tamaño del sistema") {
    assert(umbralConveniente(aleatorio(123, 1), 1.0, Nil) == 123)
    assert(umbralConveniente(Array(), 1.0, Nil) == 0)
  }
  test("P5: un solo candidato es el elegido") {
    assert(umbralConveniente(aleatorio(2000, 1), 1.0, List(500)) == 500)
  }
  test("P5: el elegido pertenece a la lista de candidatos") {
    val cand = List(100, 500, 2000)
    assert(cand.contains(umbralConveniente(aleatorio(2000, 2), 1.0, cand)))
  }
}