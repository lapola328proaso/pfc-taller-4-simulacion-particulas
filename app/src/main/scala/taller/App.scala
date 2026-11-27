package taller

import org.scalameter._

object App {

  val sim = new SimulacionParticulas()
  import sim._

  val dt = 0.5

  // Sistema determinista, sin mutación ni ciclos.
  def generar(n: Int): Sistema =
    Array.tabulate(n) { i =>
      Particula(
        Vector2(i * 0.5, (i % 100).toDouble),
        Vector2((i % 7) - 3.0, (i % 5) - 2.0),
        1.0 + (i % 3)
      )
    }

  // Mide dos bloques ALTERNADOS (a, b, a, b...) para que el ruido de la máquina les
  // afecte por igual. Libera memoria antes de cada corrida, descarta las 4 primeras
  // (calentamiento) y devuelve la mediana en ms de cada bloque.
  def medirPar(a: => Any, b: => Any): (Double, Double) = {
    def tiempo(f: => Any): Double = {
      System.gc()
      val t0 = System.nanoTime()
      val _ = f
      (System.nanoTime() - t0) / 1e6
    }
    def mediana(xs: List[Double]): Double = xs.sorted.apply(xs.length / 2)
    val corridas = List.fill(21)((tiempo(a), tiempo(b))).drop(4)
    (mediana(corridas.map(_._1)), mediana(corridas.map(_._2)))
  }

  // Amdahl: S = 1 / (fs + (1 - fs)/p)  =>  fs = (p/S - 1) / (p - 1)
  def fraccionSerial(aceleracionMedida: Double, p: Int): Double =
    if (p <= 1) 1.0 else (p / aceleracionMedida - 1.0) / (p - 1.0)

  def main(args: Array[String]): Unit = {
    val p = Runtime.getRuntime.availableProcessors
    val tamanos = List(1000000, 2000000, 3000000) // si hay OutOfMemoryError, baje el último
    println(s"Núcleos disponibles: $p")
    println("| n | umbral | t secuencial (ms) | t paralelo (ms) | aceleración medida | aceleración estimada | fracción serial (Amdahl) |")
    println("|---|---|---|---|---|---|---|")

    tamanos.foreach { n =>
      val s = generar(n)
      val umbrales = List(100, 1000, 10000, 100000, n) // el último, n, es la fila de control
      umbrales.foreach { u =>
        val (tSec, tPar) = medirPar(pasoSecuencial(s, dt), pasoParalelo(s, dt, u))
        val acel = aceleracion(tSec, tPar)
        val est = n.toDouble / tiempoEstimado(n, p, u)
        val etiqueta = if (u == n) s"$u (control)" else u.toString
        println(f"| $n | $etiqueta | $tSec%.2f | $tPar%.2f | $acel%.2f | $est%.2f | ${fraccionSerial(acel, p)}%.3f |")
      }
    }
  }
}