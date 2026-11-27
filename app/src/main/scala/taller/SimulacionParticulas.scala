package taller

import common._
import org.scalameter._

/**
 * Taller 4: simulación de partículas repartida entre varios núcleos.
 */
class SimulacionParticulas {

  type Sistema = Array[Particula]

  // ---------------------------------------------------------------------------
  // Punto 1: un paso de la simulación
  // ---------------------------------------------------------------------------

  /** Avanza una partícula un paso dt: la posición cambia, la velocidad y la masa no. */
  def mover(p: Particula, dt: Double): Particula =
    p.copy(pos = p.pos + p.vel * dt)

  /** Calcula solo el tramo [desde, hasta) del sistema. Devuelve un arreglo de hasta - desde elementos. */
  def pasoSegmento(s: Sistema, dt: Double, desde: Int, hasta: Int): Sistema = {
    val d = math.max(0, desde)
    val h = math.min(s.length, hasta)
    Array.tabulate(math.max(0, h - d))(i => mover(s(d + i), dt))
  }

  /** Mismo primitivo que el paralelo: el sistema completo es el segmento [0, n). */
  def pasoSecuencial(s: Sistema, dt: Double): Sistema =
    pasoSegmento(s, dt, 0, s.length)

  def pasoParalelo(s: Sistema, dt: Double, umbral: Int): Sistema = {
    val u = math.max(1, umbral) // un umbral < 1 nunca terminaría de partir
    def rec(desde: Int, hasta: Int): Sistema =
      if (hasta - desde <= u) pasoSegmento(s, dt, desde, hasta)
      else {
        val medio = desde + (hasta - desde) / 2
        val (izq, der) = parallel(rec(desde, medio), rec(medio, hasta))
        izq ++ der
      }
    rec(0, s.length)
  }

  // ---------------------------------------------------------------------------
  // Punto 2: magnitudes globales (secuencial)
  // ---------------------------------------------------------------------------

  def energiaSecuencial(s: Sistema): Double = {
    def energia(p: Particula): Double = {
      val v = p.vel.norma
      0.5 * p.masa * v * v
    }
    s.foldLeft(0.0)((acc, p) => acc + energia(p))
  }

  def masaTotal(s: Sistema): Double =
    s.foldLeft(0.0)((acc, p) => acc + p.masa)

  def centroDeMasaSecuencial(s: Sistema): Vector2 = {
    val (ponderada, masa) =
      s.foldLeft((Vector2(0, 0), 0.0)) { case ((suma, m), p) =>
        (suma + p.pos * p.masa, m + p.masa)
      }
    if (masa == 0.0) Vector2(0, 0)
    else Vector2(ponderada.x / masa, ponderada.y / masa)
  }

  // ---------------------------------------------------------------------------
  // Punto 3: magnitudes globales (paralelo)
  // ---------------------------------------------------------------------------

  def energiaParalela(s: Sistema, umbral: Int): Double = {
    val u = math.max(1, umbral)
    def energia(p: Particula): Double = {
      val v = p.vel.norma
      0.5 * p.masa * v * v
    }
    def segmento(desde: Int, hasta: Int): Double =
      (desde until hasta).foldLeft(0.0)((acc, i) => acc + energia(s(i)))
    def rec(desde: Int, hasta: Int): Double =
      if (hasta - desde <= u) segmento(desde, hasta)
      else {
        val medio = desde + (hasta - desde) / 2
        val (a, b) = parallel(rec(desde, medio), rec(medio, hasta))
        a + b
      }
    rec(0, s.length)
  }

  def centroDeMasaParalelo(s: Sistema, umbral: Int): Vector2 = {
    val u = math.max(1, umbral)
    // Se reparte el par (suma ponderada de posiciones, masa), NO el promedio.
    def segmento(desde: Int, hasta: Int): (Vector2, Double) =
      (desde until hasta).foldLeft((Vector2(0, 0), 0.0)) { case ((suma, m), i) =>
        (suma + s(i).pos * s(i).masa, m + s(i).masa)
      }
    def rec(desde: Int, hasta: Int): (Vector2, Double) =
      if (hasta - desde <= u) segmento(desde, hasta)
      else {
        val medio = desde + (hasta - desde) / 2
        val ((sa, ma), (sb, mb)) = parallel(rec(desde, medio), rec(medio, hasta))
        (sa + sb, ma + mb)
      }
    val (ponderada, masa) = rec(0, s.length)
    // Una sola división, al final.
    if (masa == 0.0) Vector2(0, 0)
    else Vector2(ponderada.x / masa, ponderada.y / masa)
  }

  // ---------------------------------------------------------------------------
  // Punto 4: trabajo y profundidad
  // ---------------------------------------------------------------------------

  def trabajo(n: Int): Int = n

  /** Veces que se parte a la mitad hasta quedar en <= umbral, más el trabajo del tramo que queda.
   * Se sigue la mitad más grande (ceil), que es el camino crítico y coincide con pasoParalelo. */
  def profundidad(n: Int, umbral: Int): Int = {
    val u = math.max(1, umbral)
    def rec(m: Int): Int =
      if (m <= u) m
      else 1 + rec((m + 1) / 2)
    rec(n)
  }

  def tiempoEstimado(n: Int, p: Int, umbral: Int): Double =
    profundidad(n, umbral) + trabajo(n).toDouble / math.max(1, p)

  // ---------------------------------------------------------------------------
  // Punto 5: medir la aceleración
  // ---------------------------------------------------------------------------

  def aceleracion(tSecuencial: Double, tParalelo: Double): Double =
    tSecuencial / tParalelo

  def umbralConveniente(s: Sistema, dt: Double, umbrales: List[Int]): Int = {
    def tiempo(u: Int): Double =
      (withWarmer(new Warmer.Default) measure { pasoParalelo(s, dt, u) }).value
    if (umbrales.isEmpty) s.length
    else umbrales.map(u => (u, tiempo(u))).minBy(_._2)._1
  }
}

