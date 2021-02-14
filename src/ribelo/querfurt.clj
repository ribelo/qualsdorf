(ns ribelo.querfurt
  (:require
   [uncomplicate.fluokitten.core :as fk]
   [uncomplicate.fluokitten.jvm]
   [ribelo.halle :as h]
   [ribelo.kemnath :as math]
   [ribelo.stade :as st]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(comment
  (do (require '[criterium.core :refer [quick-bench]])
      (def data (vec (repeatedly 100000 #(/ (- 0.5 ^double (rand)) 10.0))))
      (def arr  (double-array data))))

(defn ann-return-geometric
  ^double [^long freq ret]
  (let [^doubles axs    (h/seq->double-array ret)
        n      (alength ^doubles axs)
        return (fk/fold * (fk/fmap #(+ ^double % 1.0) axs))]
    (- (math/pow return (/ freq n)) 1.0)))

(defn ann-return-simple
  ^double [^long freq ret]
  (let [^doubles arr (h/seq->double-array ret)
        mean (st/mean arr)]
    (* mean freq)))

(defn annualized-return
  "Average annualized returns over a period, convenient when comparing returns.
  It can be an Arithmetic or Geometric (default) average return: if compounded with itself the
  geometric average will be equal to the cumulative return"
  ^double [^long freq mode ret]
  (case mode
    :geometric (ann-return-geometric freq ret)
    :simple    (ann-return-simple freq ret)))

;; (defn active-return ;; TODO
;;   "Asset/Portfolio annualized return minus Benchmark annualized return"
;;   ([xs freq mode]
;;    (let [xs' (sequence (annualized-return freq mode) xs)]
;;      (comp
;;       (annualized-return freq mode)
;;       (emath/sub xs')))))

(defn annualized-risk
  "Annualized standard deviation of asset/portfolio returns"
  ^double [^long freq ret]
  (let [arr (h/seq->double-array ret)]
    (* (st/std ^doubles arr) (math/sqrt freq))))

(defn sharpe-ratio
  "Sharpe Ratio.Compute Sharpe ratio for an collection XS of values (daily, weekly, etc) and
   a free-risk rate. Annual free-risk must be divided to match the right timeframe."
  ^double [^double frisk ^long freq ret]
  (let [^doubles arr (-> (h/seq->double-array ret) (h/take-last freq))
        mean         (st/mean arr)
        std          (st/std arr)]
    (/ (- mean frisk) std)))

(defn annualized-sharpe-ratio
  ^double [^double frisk ^long freq ret]
  (let [ann-ret (ann-return-geometric freq ret)
        std     (annualized-risk freq ret)]
    (/ (- ann-ret frisk) std)))

;; (defn adjusted-sharpe-ratio ;;TODO
;;   "Sharpe Ratio adjusted for skewness and kurtosis with a penalty factor
;;    for negative skewness and excess kurtosis."
;;   ([^double frisk]
;;    (comp
;;     (x/transjuxt [(sharpe-ratio frisk) st/skewness st/kurtosis])
;;     (x/reduce
;;      (fn
;;        ([] (transient []))
;;        ([[sr sk ku]] (* sr (- (+ 1 (* (/ sk 6) sr))
;;                               (* (/ (- ku 3) 24) (math/sqrt sr)))))
;;        ([acc [sr sk ku]] (-> acc (conj! sr) (conj! sk) (conj! ku)))))))
;;   ([] (adjusted-sharpe-ratio 0.0)))

;; (defn annualized-adjusted-sharpe-ratio ;;TODO
;;   "Sharpe Ratio adjusted for skewness and kurtosis with a penalty factor
;;    for negative skewness and excess kurtosis."
;;   ([^double frisk ^long freq mode]
;;    (comp
;;     (x/transjuxt [st/skewness st/kurtosis
;;                   (annualized-return freq mode) (annualized-risk freq)])
;;     (x/reduce
;;      (fn
;;        ([] (transient []))
;;        ([[sk ku annret annrisk]]
;;         (let [sr (/ (- (/ (Math/round ^double (* 10000 annret)) 10000) frisk)
;;                     (/ (Math/round ^double (* 10000 annrisk)) 10000))]
;;           (* sr (- (+ 1 (* (/ sk 6) sr))
;;                    (* (- ku 3) (math/sqrt sr))))))
;;        ([acc [sk ku annret annrisk]]
;;         (-> acc (conj! sk) (conj! ku) (conj! annret) (conj! annrisk)))))))
;;   ([^double frisk ^long freq]
;;    (annualized-adjusted-sharpe-ratio frisk freq :geometric))
;;   ([^double frisk]
;;    (annualized-adjusted-sharpe-ratio frisk 252 :geometric))
;;   ([]
;;    (annualized-adjusted-sharpe-ratio 0.0 252 :geometric)))

(defn downside-risk
  "Downside Risk or Semi-Standard Deviation.
   Measures the variability of underperformance below a minimum target rate"
  ^double [^double mar ret]
  (let [^doubles arr (h/seq->double-array ret)
        n   (alength ^doubles arr)]
    (loop [i 0 sum 0.0]
      (if (< i n)
        (recur (inc i)
               (+ sum
                  (/ (math/sq (math/min 0.0 (- (aget ^doubles arr i) mar))) n)))
        (math/sqrt sum)))))

(defn sortino-ratio
  "Sortino ratio"
  ^double [^double frisk ^double mar ret]
  (let [^doubles arr (h/seq->double-array ret)
        dr           (downside-risk mar arr)
        mean         (st/mean arr)]
    (/ (- mean frisk) dr)))

(defn drawdown
  "Drawdowon from Peak. Any continuous losing return period."
  ^doubles [ret]
  (let [^doubles arr (h/seq->double-array ret)
        n            (alength ^doubles arr)
        r            (double-array n)]
    (loop [i 0 s 1.0 mx 1.0]
      (if (< i n)
        (let [v   (* s (+ 1.0 (aget ^doubles arr i)))
              mx' (math/max v mx)
              dr  (/ (- mx' v) mx')]
          (aset r i dr)
          (recur (inc i) v mx'))
        r))))

(defn continuous-drawdown
  ^doubles [ret]
  (let [^doubles arr (h/seq->double-array ret)
        n            (alength arr)
        dq           (java.util.ArrayDeque.)]
    (loop [i 0 s 1.0]
      (when (< i n)
        (let [v (aget ^doubles arr i)]
          (cond
            (and (zero? i) (< v 0.0))
            (recur (inc i) (+ 1.0 v))
            (and (zero? i) (> v 0.0))
            (recur (inc i) 1.0)
            (< 0 i)
            (if (< v 0.0)
              (recur (inc i) (* s (+ 1.0 v)))
              (let [dd (- 1.0 s)]
                (when-not (zero? dd)
                  (.add dq dd))
                (recur (inc i) 1.0)))))))
    (let [r (double-array (.toArray dq))]
      (.clear dq)
      r)))

(comment
  (quick-bench (continuous-drawdown data)))

(defn average-drawdown
  ^double [ret]
  (let [^doubles arr (h/seq->double-array ret)]
    (->> (continuous-drawdown arr)
         (st/mean))))

(defn maximum-drawdown
  ^double [ret]
  (let [^doubles arr (h/seq->double-array ret)]
    (->> (continuous-drawdown arr)
         (st/max))))

(defn rate-of-return
  "Simple rate of return calculated from the last and the first value of
  an array of numbers."
  ^double [ret]
  (let [^doubles arr (as-> (h/seq->double-array ret) $
                       (fk/fmap #(+ 1.0 ^double %) $)
                       (h/reductions $ *))
        l            (h/last arr)
        f            (h/first arr)]
    (- (/ ^double l ^double f) 1.0)))

(defn rate-of-change
  "Simple rate of chane calculated from the last and the first value of
  an array of numbers."
  ^double [^long n ret]
  (let [^doubles arr (h/seq->double-array ret)
        c            (alength arr)]
    (if (< n c)
      (let [l (h/last arr)
            i (aget ^doubles arr (dec (- c n)))]
        (- (/ ^double l i) 1.0))
      0.0)))

(defn cagr
  "Compound annual growth rate"
  ^double [^double n ret]
  (let [^doubles arr (h/seq->double-array ret)]
    (- (math/pow
        (+ 1.0
           (rate-of-return arr))
        (/ 1.0 n))
       1.0)))

(defn calmar-ratio
  "A risk-adjusted measure like Sharpe ratio that uses maximum drawdown instead of
  standard deviation for risk."
  [^double frisk ^long freq ret]
  (let [^doubles arr (h/seq->double-array ret)
        maxdd        (maximum-drawdown ret)
        annret       (ann-return-geometric freq arr)]
    (/ (- annret frisk) maxdd)))

;; (defn downside-potential ;;TODO
;;   "Downside potential"
;;   ([mar]
;;    (comp (x/transjuxt [(comp (emath/mul -1.0) (emath/add mar) (emath/max 0.0) (x/into [])) x/count])
;;          (x/reduce
;;           (fn
;;             ([] (transient []))
;;             ([[xs count]] (transduce (emath/div count) + xs))
;;             ([acc [xs count]] (-> acc (conj! xs) (conj! count)))))))
;;   ([] (downside-potential 0.0)))

;; (defn burke-ratio ;;TODO
;;   "A risk-adjusted measure with free risk and drawdowns.
;;    For the 'simple' mode the excess return over free risk is divided by the square root of
;;    the sum of the square of the drawdowns. For the 'modified' mode the Burke Ratio is multiplied
;;    by the square root of the number of datas."
;;   ([frisk freq mode]
;;    (comp
;;     (x/transjuxt [(annualized-return freq)
;;                   (comp continuous-drawdown
;;                         (emath/pow 2)
;;                         (x/reduce +)
;;                         emath/sqrt)
;;                   x/count])
;;     (x/reduce
;;      (fn
;;        ([] (transient []))
;;        ([[annret dd c]] (case mode
;;                           :simple (/ (- annret frisk) dd)
;;                           :modified (* (/ (- annret frisk) dd)
;;                                        (math/sqrt c))))
;;        ([acc [annret dd c]] (-> acc
;;                                 (conj! annret)
;;                                 (conj! dd)
;;                                 (conj! c)))))))
;;   ([frisk freq]
;;    (burke-ratio frisk freq :simple))
;;   ([frisk]
;;    (burke-ratio frisk 252 :simple))
;;   ([]
;;    (burke-ratio 0.0 252 :simple)))

;; (def ulcer-index ;;TODO
;;   "Ulcer Index of Peter G. Martin (1987). The impact of long, deep drawdowns will have significant
;;   impact because the underperformance since the last peak is squared."
;;   (comp
;;    (x/transjuxt [(comp drawdown
;;                        (emath/pow 2)
;;                        (x/reduce +))
;;                  x/count])
;;    (x/reduce
;;     (fn
;;       ([] (transient []))
;;       ([[dd c]] (math/sqrt (/ dd c)))
;;       ([acc [dd c]] (-> acc (conj! dd) (conj! c)))))))

;; (defn martin-ratio ;;TODO
;;   "A risk-adjusted measure with free risk and Ulcer index.
;;    Martin Ratio = (Portfolio Return - RiskFree) / Ulcer Index
;;    Mode: :return, :geometric (default: :return)"
;;   ([frisk freq]
;;    (comp
;;     (x/transjuxt [(annualized-return freq)
;;                   ulcer-index])
;;     (x/reduce
;;      (fn
;;        ([] (transient []))
;;        ([[annret u]] (/ (- annret frisk) u))
;;        ([acc [annret u]] (-> acc (conj! annret) (conj! u)))))))
;;   ([frisk]
;;    (martin-ratio frisk 252))
;;   ([]
;;    (martin-ratio 0.0 252)))

;; (def hurst-index ;;TODO
;;   "It's a useful statistic for detecting if a time series is mean reverting (anti-persistent), totally random or persistent.
;;    A value in the range [0.5) indicates mean-reverting (anti-persistent)
;;    A value of 0.5 indicate a random walk
;;    A value H in the range (0.5,1] indicates momentum (persistent)"
;;   (comp
;;    (x/transjuxt [(comp emath/cumdev st/max)
;;                  (comp emath/cumdev st/min)
;;                  st/std
;;                  x/count])
;;    (x/reduce
;;     (fn
;;       ([] (transient []))
;;       ([[mx mn std n]]
;;        (let [rs (/ (- mx mn) std)]
;;          (/ (math/log rs) (math/log n))))
;;       ([acc [mx mn std n]] (-> acc
;;                                (conj! mx) (conj! mn)
;;                                (conj! std) (conj! n)))))))

;; (defn info-ratio ;;TODO
;;   "Information Ratio"
;;   [xs]
;;   (comp
;;    (x/transjuxt [(comp (emath/sub xs) st/std)
;;                  (comp (emath/sub xs) st/mean)])
;;    (x/reduce
;;     (fn
;;       ([] (transient []))
;;       ([[std mean]] (/ mean std))
;;       ([acc [std mean]] [std mean])))))

;; (defn jensen-alpha ;;TODO
;;   "Ex-post alpha calculated with regression line.
;;   Free-risk is the avereage free-risk for the timeframe selected."
;;   [xs frisk]
;;   (comp
;;    (x/transjuxt [(comp (emath/sub xs) st/std)
;;                  (comp (emath/sub xs) st/mean)])
;;    (x/reduce
;;     (fn
;;       ([] (transient []))
;;       ([[std mean]]
;;        (/ mean std))
;;       ([acc [std mean]] (-> acc (conj! std) (conj! mean)))))))

;; (defn modigliani ;;TODO
;;   "Modigliani index for risk-adjusted return"
;;   [ys frisk]
;;   (let [[stdb] (sequence st/std ys)]
;;     (comp
;;      (x/transjuxt [st/mean
;;                    (sharpe-ratio frisk)
;;                    st/std])
;;      (x/reduce
;;       (fn
;;         ([] [0.0 0.0 0.0])
;;         ([[mean sharpe std]]
;;          (+ mean (* sharpe (- stdb std))))
;;         ([acc coll] coll))))))

(defn rolling-economic-drawndown
  ^double [^long freq ret]
  (let [^doubles arr (-> (h/seq->double-array ret) (h/take-last freq))
        mx           (st/max ^doubles arr)]
    (- 1.0 (/ ^double (h/last arr) mx))))

(defn tick->ret
  "Convert a value series to a return series"
  ^doubles [close]
  (let [^doubles arr (h/seq->double-array close)
        n            (alength arr)
        nn           (dec n)
        r            (double-array (dec n))]
    (loop [i 0]
      (when (< i nn)
        (let [x   (aget ^doubles arr i)
              y   (aget ^doubles arr (inc i))
              ret (- (/ y x) 1.0)]
          (aset ^doubles r i ret)
          (recur (inc i)))))
    r))

(defn redp-single-allocation
  ^double [^double frisk ^double risk ^long freq close]
  (let [^doubles close* (h/seq->double-array close)
        ret             (tick->ret close*)
        redp            (rolling-economic-drawndown freq close*)
        std             (annualized-risk freq ret)
        sr              (annualized-sharpe-ratio frisk freq ret)]
    (math/min 1.0
              (math/max 0.0 (* (/ (+ (/ sr std) 0.5)
                                  (- 1.0 (math/pow risk 2.0)))
                               (math/max 0.0 (/ (- risk redp)
                                                (- 1.0 redp))))))))

(defn- redp-stats [^double frisk ^double risk ^long freq ^doubles close]
  (let [close*  (h/seq->double-array close)
        ret     (tick->ret close*)
        std     (annualized-risk freq ret)
        ann-ret (ann-return-geometric freq ret)
        drift   (math/max 0.0 (+ (- ann-ret frisk) (/ (math/sq std) 2.0)))
        redp    (rolling-economic-drawndown freq close*)
        Y       (* (/ 1.0 (- 1.0 (math/sq risk)))
                   (/ (- risk redp)
                      (- 1.0 redp)))]
    [std drift Y]))

(defn redp-multiple-allocation
  ([frisk risk freq assets]
   (let [[volatility
          drift
          Y]            (apply map vector (fk/fmap (partial redp-stats frisk risk freq) assets))
         inv-volatility (fk/fmap #(math/pow % -1.0) volatility)
         matrix         (fk/fmap (comp math/round *) inv-volatility drift inv-volatility Y)
         sum            (fk/fold + matrix)]
     (if (<= ^double sum 1.0) matrix (fk/fmap #(/ ^double % ^double sum) matrix)))))
