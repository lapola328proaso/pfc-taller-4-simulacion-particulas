# Informe de corrección

Notación: el sistema es $s = (s_0, \dots, s_{n-1})$ con $s_i = (x_i, v_i, m_i)$. Un tramo es el intervalo $[a, b)$ de índices. Para $a<b$ se define $\mu(a,b) = a + \lfloor (b-a)/2 \rfloor$, el `medio` del código.

## 1. Punto 1: un paso de la simulación

**Especificación.** $\text{paso}(s, dt)_i = (x_i + v_i\,dt,\ v_i,\ m_i)$ para todo $0 \le i < n$.

**`pasoSegmento(s, dt, a, b)`** devuelve el arreglo $r$ de longitud $b-a$ con $r_j = \text{paso}(s,dt)_{a+j}$, por construcción con `Array.tabulate`. Por tanto `pasoSecuencial` $=$ `pasoSegmento(s, dt, 0, n)` cumple la especificación.

### Lema 1: partir en segmentos no pierde ni repite partículas

Sea $a<b$ y $\mu = \mu(a,b)$. Entonces $a < \mu < b$ cuando $b-a \ge 2$, y

$$[a,b) = [a,\mu) \,\cup\, [\mu,b), \qquad [a,\mu) \cap [\mu,b) = \varnothing .$$

*Prueba.* Todo $i$ con $a \le i < \mu$ no cumple $i \ge \mu$ y viceversa, así que los intervalos son disjuntos. Todo $i$ en $[a,b)$ cumple $i<\mu$ o $i\ge\mu$, así que la unión es el intervalo completo. $\mu>a$ porque $\lfloor (b-a)/2\rfloor \ge 1$ si $b-a\ge2$, y $\mu<b$ porque $\lfloor (b-a)/2\rfloor < b-a$. Por eso ambas llamadas recursivas reciben un tramo estrictamente más corto. $\square$

### Teorema 1: `pasoParalelo` $=$ `pasoSecuencial`

Se prueba por inducción fuerte sobre $\ell = b-a$ que `rec(a,b)` $=$ `pasoSegmento(s,dt,a,b)`, para cualquier umbral $u\ge1$ (el código usa $\max(1,u)$).

- **Base** ($\ell \le u$): `rec` llama directamente a `pasoSegmento`.
- **Paso** ($\ell > u \ge 1$, luego $\ell\ge2$): por el Lema 1, los tramos $[a,\mu)$ y $[\mu,b)$ son más cortos. Por hipótesis, `rec(a,μ)` y `rec(μ,b)` valen `pasoSegmento` de esos tramos. Como cada uno calcula el mismo valor ya sea en secuencial o en otra tarea (las funciones solo leen $s$ y no hay estado compartido modificado), el par que devuelve `parallel` es el mismo. La concatenación $\text{izq}\,\Vert\,\text{der}$ tiene en las posiciones $0..\mu-a-1$ los elementos $\text{paso}_{a}, \dots, \text{paso}_{\mu-1}$ y a continuación $\text{paso}_{\mu},\dots,\text{paso}_{b-1}$: exactamente `pasoSegmento(s,dt,a,b)`.

Con $a=0$, $b=n$ se obtiene el resultado. Vale para todo $u$ y todo $dt$. Con $u\ge n$ no hay llamada a `parallel`, y el cálculo coincide con `pasoSecuencial` operación por operación. $\blacksquare$

## 2. Punto 2: magnitudes globales

- $E(s) = \sum_{i} \tfrac12 m_i \lVert v_i \rVert^2$. `energiaSecuencial` es el `foldLeft` de esa suma; el vacío da $0$.
- $M(s) = \sum_i m_i$.
- $C(s) = \dfrac{\sum_i m_i\,x_i}{M(s)}$ si $M(s)\ne0$, y $C(s)=(0,0)$ si $M(s)=0$. Se acumula el par $(S,M)=(\sum m_ix_i,\ \sum m_i)$ y se divide una vez.

### Corrección del recorrido secuencial de un tramo (invariante)

Las hojas de `energiaParalela` y `centroDeMasaParalelo` recorren el tramo $[a,b)$ con un `foldLeft`. Se argumenta como un proceso iterativo:

- Estado: $s=(i,\ acc)$, con $acc$ la energía acumulada (para el centro de masa, el par $(S,M)$).
- Estado inicial: $s_0=(a,\ 0)$.
- Estado final: $i=b$.
- Transformación: $\text{transformar}(i,acc)=(i+1,\ acc+e(s_i))$, con $e(s_i)=\tfrac12 m_i\lVert v_i\rVert^2$.
- Invariante: $\text{Inv}(i,acc)\equiv a\le i\le b\ \land\ acc=\sum_{j=a}^{i-1}e(s_j)$.

**1.** $\text{Inv}(s_0)$: con $i=a$ el rango $j=a,\dots,a-1$ es vacío, y $acc=0$ es la suma vacía.

**2.** Si $\text{Inv}(i,acc)$ e $i<b$, entonces el nuevo estado es $(i+1,\ acc+e(s_i))$. Se cumple $a\le i+1\le b$ y
$$acc+e(s_i)=\sum_{j=a}^{i-1}e(s_j)+e(s_i)=\sum_{j=a}^{i}e(s_j).$$

**3.** En el estado final $i=b$, la invariante da $acc=\sum_{j=a}^{b-1}e(s_j)$, que es lo que el tramo debe devolver.

**4.** Terminación: en cada paso $b-i$ baja en 1, así que se llega a $i=b$ tras $b-a$ pasos.

El mismo argumento vale para el centro de masa cambiando $e(s_i)$ por el par $(m_i x_i,\ m_i)$ y la suma por $\oplus$.

## 3. Punto 3: por qué la suma sí se reparte y el promedio no

### 3.1 La operación de combinación es asociativa (y tiene neutro)

Para `energiaParalela` se combina con $+$ en $\mathbb{R}$. Para `centroDeMasaParalelo` se combina con $\oplus$ sobre pares $(S,M)\in\mathbb{R}^2\times\mathbb{R}$:

$$(S_1,M_1)\oplus(S_2,M_2) = (S_1+S_2,\ M_1+M_2).$$

Es suma componente a componente, y la suma de reales es asociativa y conmutativa:

$$\big((S_1,M_1)\oplus(S_2,M_2)\big)\oplus(S_3,M_3) = (S_1+S_2+S_3,\ M_1+M_2+M_3) = (S_1,M_1)\oplus\big((S_2,M_2)\oplus(S_3,M_3)\big).$$

El neutro es $((0,0),0)$ (tramo vacío). Por eso no importa cómo se agrupen los tramos: cualquier árbol de partición da la misma suma total, que es la del secuencial.

### 3.2 El promedio no es combinable directamente

Sean dos tramos con centros $c_1,c_2$ y masas $M_1,M_2$. El centro del conjunto es

$$C = \frac{S_1+S_2}{M_1+M_2} = \frac{M_1 c_1 + M_2 c_2}{M_1+M_2},$$

que solo es el promedio simple $\frac{c_1+c_2}{2}$ cuando $M_1=M_2$. Además, el "promedio de dos promedios" no es asociativo: $\operatorname{prom}(\operatorname{prom}(a,b),c)\ne\operatorname{prom}(a,\operatorname{prom}(b,c))$ en general: con $a=0,b=0,c=12$ el primero da $6$ y el segundo da $3$.

**Contraejemplo con partículas.** $P_1$ en $(0,0)$ con masa $1$ y $P_2$ en $(10,0)$ con masa $3$. Promediando cada mitad: $(0+10)/2=5$. Con la suma ponderada: $(0\cdot1+10\cdot3)/(1+3)=7.5$. El valor correcto es $7.5$; el promedio de promedios da $5$. Por eso se reparte el par $(S,M)$ y se divide una sola vez en la raíz.

### 3.3 Teorema 2: las versiones paralelas coinciden con las secuenciales

Por inducción sobre $\ell=b-a$, igual que el Teorema 1: `rec(a,b)` devuelve $\sum_{i=a}^{b-1} f(s_i)$ con $f$ la contribución de cada partícula (energía, o el par $(m_ix_i,m_i)$).

- Base ($\ell\le u$): el segmento recorre justo $a,\dots,b-1$.
- Paso: por el Lema 1, $\sum_{i=a}^{b-1} f(s_i) = \sum_{i=a}^{\mu-1} f(s_i) + \sum_{i=\mu}^{b-1} f(s_i)$, y por hipótesis cada sumando es lo que devuelve la rama correspondiente. La combinación $+$ o $\oplus$ es asociativa (3.1), así que el resultado es la suma sobre $[0,n)$. En el centro de masa se divide después, una vez, entre $M$.

### 3.4 Aritmética de punto flotante

La asociatividad vale en $\mathbb{R}$, pero no exactamente en los `Double` de IEEE 754: sumar en otro orden puede cambiar los últimos dígitos. Por eso las pruebas comparan las versiones paralelas con una **tolerancia** y no con igualdad exacta. En `pasoParalelo` no hay este problema: no se suma nada, cada posición se calcula con las mismas operaciones, y la igualdad es exacta.

## 4. Punto 4: trabajo y profundidad

- $W(n)=n$: una operación por partícula. Partir no cambia cuántas se hacen, solo dónde.
- $D(n,u)=n$ si $n\le u$, y $D(n,u)=1+D(\lceil n/2\rceil,u)$ si $n>u$. Termina porque $\lceil n/2\rceil<n$ para $n\ge2$ (y $u\ge1$ garantiza $n\ge2$ al partir). Sigue la mitad mayor, que es la que tarda más al correr en paralelo.
- Para $n$ potencia de 2 y $u$ también, hay $\log_2(n/u)$ particiones y queda un tramo de $u$: $D=\log_2(n/u)+u$. Ejemplo: $D(1024,1)=10+1=11$ y $D(1024,512)=1+512=513$.
- $T_p(n,u)=D(n,u)+W(n)/p$.

## 5. Punto 5

- `aceleracion` es el cociente $S=T_{sec}/T_{par}$; $S=1$ si ambos tiempos son iguales.
- `umbralConveniente` devuelve el $u^\*=\arg\min_{u\in U} t(u)$ sobre los candidatos $U$ y $|s|$ si $U=\varnothing$. Como se mide con calentamiento de la JVM, el resultado es el de menor tiempo medido; pequeñas diferencias entre umbrales cercanos pueden ser ruido de la máquina.