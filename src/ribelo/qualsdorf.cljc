(ns ribelo.qualsdorf
  (:require
   [ribelo.halle :as h]
   [ribelo.kemnath :as math]
   [ribelo.stade :as stats]))

#?(:clj (set! *warn-on-reflection* true))
#?(:clj (set! *unchecked-math* :warn-on-boxed))

(comment
  (do (require '[taoensso.encore :as enc])
      (def data (vec (repeatedly 100000 #(/ (- 0.5 ^double (rand)) 10.0))))
      (def arr  (double-array data))
      (def X [0.003,0.026,0.015,-0.009,0.014,0.024,0.015,0.066,-0.014,0.039])))

(defn annualized-return-geometric
  (^double [           ret] (annualized-return-geometric 252 ret))
  (^double [^long freq ret]
   (let [arr    (h/->double-array ret)
         n      (alength arr)
         return (h/reduce * (h/map #(+ ^double % 1.0) arr))]
     (- (math/pow return (/ freq n)) 1.0))))

(comment
  (annualized-return-geometric 12 X))

(defn ann-return-simple
  (^double [           ret] (ann-return-simple 252 ret))
  (^double [^long freq ret]
   (let [arr  (h/->double-array ret)
         mean (stats/mean arr)]
     (* mean freq))))

(defn annualized-return
  "average annualized returns over a period, convenient when comparing returns.
  it can be an arithmetic or geometric (default) average return: if compounded with itself the
  geometric average will be equal to the cumulative return"
  (^double [                ret] (annualized-return 252  :geometric ret))
  (^double [^long freq      ret] (annualized-return freq :geometric ret))
  (^double [^long freq mode ret]
   (case mode
     :geometric (annualized-return-geometric freq ret)
     :simple    (ann-return-simple freq ret))))

(defn active-return
  "asset/portfolio annualized return minus benchmark annualized return"
  (^double [                r1 r2] (active-return 252  :geometric r1 r2))
  (^double [^long freq      r1 r2] (active-return freq :geometric r1 r2))
  (^double [^long freq mode r1 r2]
   (- (annualized-return freq mode r1)
      (annualized-return freq mode r2))))

(defn excess-return
  "return on asset - risk free rate"
  (^doubles [              ret] (excess-return 0.0 ret))
  (^doubles [^double frisk ret]
   (let [arr (h/->double-array ret)]
     (stats/sub arr frisk))))

(defn annualized-risk
  "annualized standard deviation of asset/portfolio returns"
  ^double [^long freq ret]
  (let [arr (->> (h/->double-array ret) (h/take-last freq))]
    (* (stats/std arr) (math/sqrt freq))))

(defn sharpe-ratio
  "sharpe ratio.compute sharpe ratio for an collection xs of values (daily, weekly, etc) and
   a free-risk rate. annual free-risk must be divided to match the right timeframe."
  (^double [              ret] (sharpe-ratio 0.0 ret))
  (^double [^double frisk ret]
   (let [arr  (h/->double-array ret)
         mean (stats/mean arr)
         std  (stats/std arr)]
     (/ (- mean frisk) std))))

(defn annualized-sharpe-ratio
  (^double [                         ret] (annualized-sharpe-ratio 0.0   (count ret) ret))
  (^double [^double frisk            ret] (annualized-sharpe-ratio frisk (count ret) ret))
  (^double [^double frisk ^long freq ret]
   (let [ann-ret  (annualized-return-geometric freq (excess-return frisk ret))
         ann-risk (annualized-risk freq ret)]
     (/ ann-ret ann-risk))))

(defn annualized-adjusted-sharpe-ratio
  "sharpe ratio adjusted for skewness and kurtosis with a penalty factor for
  negative skewness and excess kurtosis."
  (^double [                              ret]
   (annualized-adjusted-sharpe-ratio 0.0   (count ret) :geometric ret))
  (^double [^double frisk                 ret]
   (annualized-adjusted-sharpe-ratio frisk (count ret) :geometric ret))
  (^double [^double frisk            mode ret]
   (annualized-adjusted-sharpe-ratio frisk (count ret) mode       ret))
  (^double [^double frisk ^long freq mode ret]
   (let [arr   (h/->double-array ret)
         aret  (annualized-return freq mode ret)
         arisk (annualized-risk freq ret)
         sr    (/ (- aret frisk) arisk)
         sk    (stats/skewness arr)
         ku    (stats/kurtosis arr)]
     (* sr (+ 1.0 (- (* (/ sk 6) sr) (* (/ (- ku 3) 24) (math/sq sr))))))))

(defn downside-risk
  "downside risk or semi-standard deviation.
  measures the variability of underperformance below a minimum target rate"
  (^double [            ret] (downside-risk 0.0 ret))
  (^double [^double mar ret]
   (let [a1 (h/->double-array ret)
         n  (alength a1)]
     (loop [i 0 acc 0.0]
       (if (< i n)
         (recur (unchecked-inc-int i)
                (+ acc (/ (math/sq (math/min (- (aget a1 i) mar) 0.0)) n)))
         (math/sqrt acc))))))

(defn downside-potential
  "downside potential is the first lower partial moment"
  (^double [            ret] (downside-potential 0.0 ret))
  (^double [^double mar ret]
   (let [a1 (h/->double-array ret)
         n  (alength a1)]
     (loop [i 0 acc 0.0]
       (if (< i n)
         (recur (unchecked-inc-int i) (+ acc (/ (math/max (- mar (aget a1 i)) 0) n)))
         acc)))))

(defn upside-risk
  "measures the variability of overperformance above a minimum target rate"
  (^double [            ret] (downside-risk 0.0 ret))
  (^double [^double mar ret]
   (let [a1 (h/->double-array ret)
         n  (alength a1)]
     (loop [i 0 acc 0.0]
       (if (< i n)
         (recur (unchecked-inc-int i)
                (+ acc (/ (math/sq (math/min (- mar (aget a1 i)) 0.0)) n)))
         (math/sqrt acc))))))

(defn upside-potential
  (^double [            ret] (upside-potential 0.0 ret))
  (^double [^double mar ret]
   (let [a1 (h/->double-array ret)
         n  (alength a1)]
     (loop [i 0 acc 0.0]
       (if (< i n)
         (recur (unchecked-inc-int i) (+ acc (/ (math/max (- (aget a1 i) mar) 0) n)))
         acc)))))

(defn sortino-ratio
  "Sortino ratio"
  (^double [                          ret] (sortino-ratio 0.0 0.0 ret))
  (^double [              ^double mar ret] (sortino-ratio 0.0 mar ret))
  (^double [^double frisk ^double mar ret]
   (let [^doubles a1 (h/->double-array ret)
         dr          (downside-risk mar a1)
         mean        (stats/mean a1)]
     (/ (- mean frisk) dr))))

(defn annualized-sortino-ratio
  (^double [                       ret] (annualized-sortino-ratio 0.0 (count ret) ret))
  (^double [^double mar            ret] (annualized-sortino-ratio mar (count ret) ret))
  (^double [^double mar ^long freq ret]
   (/ (- (annualized-return-geometric freq ret) mar)
      (* (downside-risk (/ mar freq) ret) (math/sqrt freq)))))

(defn drawdown
  "Drawdowon from Peak. Any continuous losing return period."
  ^doubles [ret]
  (let [a1 (h/->double-array ret)
        n  (alength ^doubles a1)
        r  (double-array n)]
    (loop [i 0 s 1.0 mx 1.0]
      (if (< i n)
        (let [v   (* s (+ 1.0 (aget ^doubles a1 i)))
              mx' (math/max v mx)
              dr  (/ (- mx' v) mx')]
          (aset r i dr)
          (recur (inc i) v mx'))
        r))))

(defn continuous-drawdown
  ^doubles [ret]
  (let [a1 (h/->double-array ret)
        n  (alength a1)
        dq #?(:clj (java.util.ArrayDeque.) :cljs #js [])]
    (loop [i 0 s 1.0]
      (if (< i n)
        (let [v (aget ^doubles a1 i)]
          (cond
            (and (zero? i) (< v 0.0))
            (recur (inc i) (+ 1.0 v))
            ;;
            (and (zero? i) (> v 0.0))
            (recur (inc i) 1.0)
            ;;
            (> i 0)
            (if (< v 0.0)
              (recur (inc i) (* s (+ 1.0 v)))
              (let [dd (- 1.0 s)]
                (when-not (zero? dd)
                  #?(:clj (.add dq dd) :cljs (.push dq dd)))
                (recur (inc i) 1.0)))))
        (when (< s 1.0)
          (let [dd (- 1.0 s)]
            (when-not (zero? dd)
              #?(:clj (.add dq dd) :cljs (.push dq dd)))))))
    (let [r #?(:clj (double-array (.toArray dq)) :cljs dq)]
      #?(:clj (.clear dq))
      r)))

(defn average-drawdown
  ^double [ret]
  (stats/mean (continuous-drawdown ret)))

(defn maximum-drawdown
  ^double [ret]
  (stats/max (continuous-drawdown ret)))

(defn ulcer-index
  "ulcer index of peter g. martin (1987). the impact of long, deep drawdowns will
  have significant * impact because the underperformance since the last peak is
  squared."
  ^double [ret]
  (let [a1 (h/->double-array ret)
        dd (drawdown a1)]
    (math/sqrt (/ (stats/sum (stats/sq dd)) (alength a1)))))

(defn calmar-ratio
  "a risk-adjusted measure like sharpe ratio that uses maximum drawdown instead of
  standard deviation for risk."
  (^double [              ret] (calmar-ratio 0.0 ret))
  (^double [^double frisk ret]
   (let [a1     (h/->double-array ret)
         n      (alength a1)
         maxdd  (maximum-drawdown ret)
         annret (annualized-return-geometric n a1)]
     (/ (- annret frisk) maxdd))))

(defn hist-var
  "historical value at risk"
  (^double [                     ret] (hist-var 0.95 ret))
  (^double [^double p ret]
   (stats/quantile (- 1.0 p) (h/->double-array ret))))

(defn rate-of-return
  "simple rate of return calculated from the last and the first value of an array
  of numbers."
  ^double [ret]
  (let [a1 (as-> (h/->double-array ret) $
             (h/map #(+ 1.0 ^double %) $)
             (h/reductions * $))
        l  (h/last a1)
        f  (h/first a1)]
    (- (/ ^double l ^double f) 1.0)))

(defn rate-of-change
  "simple rate of change calculated from the last and the first value of an array
  of numbers."
  ^double [^long n ret]
  (let [a1 (h/->double-array ret)
        c  (alength a1)]
    (if (< n c)
      (let [l (h/last a1)
            i (aget a1 (dec (- c n)))]
        (- (/ ^double l i) 1.0))
      0.0)))

(defn cagr
  "Compound annual growth rate"
  ^double [ret]
  (let [a1 (h/->double-array ret)]
    (- (math/pow
        (+ 1.0
           (rate-of-return a1))
        (/ 1.0 (alength a1)))
       1.0)))

(defn rolling-economic-drawndown
  ^double [close]
  (let [a1 (h/->double-array close)
        mx (stats/max a1)]
    (- 1.0 (/ ^double (h/last a1) mx))))

(defn tick->ret
  "convert a value series to a return series"
  ^doubles [close]
  (let [arr (h/->double-array close)
        n   (alength arr)
        nn  (dec n)
        r   (h/double-array (dec n))]
    (loop [i 0]
      (when (< i nn)
        (let [x   (aget arr i)
              y   (aget arr (inc i))
              ret (- (/ y x) 1.0)]
          (aset r i ret)
          (recur (inc i)))))
    r))

(defn redp-single-allocation
  ^double [^double frisk ^double risk close]
  (let [close' (h/->double-array close)
        n      (alength close')
        ret    (tick->ret close')
        redp   (rolling-economic-drawndown close')
        std    (annualized-risk n ret)
        sr     (annualized-sharpe-ratio frisk n ret)]
    (math/min
     1.0
     (math/max
      0.0
      (* (/ (+ (/ sr std) 0.5)
            (- 1.0 (math/sq risk)))
         (math/max 0.0 (/ (- risk redp)
                          (- 1.0 redp))))))))

(defn redp-stats [^double frisk ^double risk ^doubles close]
  (let [ca      (h/->double-array close)
        n       (alength ca)
        ret     (tick->ret ca)
        std     (annualized-risk n ret)
        ann-ret (annualized-return-geometric n ret)
        drift   (math/max 0.0 (+ (- ann-ret frisk) (/ (math/sq std) 2.0)))
        redp    (rolling-economic-drawndown ca)
        Y       (* (/ 1.0 (- 1.0 (math/sq risk)))
                   (/ (- risk redp)
                      (- 1.0 redp)))]
    [std drift Y]))

(defn redp-multiple-allocation
  ([frisk risk assets]
   (let [[volatility
          drift
          Y]            (stats/transpose (mapv (partial redp-stats frisk risk) assets))
         inv-volatility (h/map #(math/pow % -1.0) volatility)
         matrix         (mapv * inv-volatility drift inv-volatility Y)
         sum            (h/reduce + matrix)]
     (if (<= ^double sum 1.0) matrix (h/map #(math/round2 (/ ^double % ^double sum)) matrix)))))

(comment
  (require '[taoensso.encore :as enc])

  (let [v (vec (range 1000))
        arr (double-array v)]
    (enc/qb 1e3
      (reduce + v)
      (h/reduce + arr)))

  (enc/qb (h/map (fn [x y] (+ x y))
                 (vec (range 1000))
                 (vec (range 1000))
                 ;; (vec (range 1000))
                 ))
  (h/map (comp math/round *)
         [1 2 3]
         [1 2 3]
         [1 2 3]))
