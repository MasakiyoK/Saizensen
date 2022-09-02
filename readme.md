version Aug. 2022 / by Masakiyo Kitazawa
文責：北沢正清

This readme is an instruction and tips for the codes placed at
https://github.com/MasakiyoK/Saizensen/


They are developed for a book 
"Quark Matter under Eetreme Conditions Phase Transitions in the World of Elementary Particles",
Masakiyo Kitazawa, Teiji Kunihiro, 2022, ISBN:978-4320035492 

# 概要

このreadmeは、  
「物理学最前線29：超高温・高密度のクォーク物質〜素粒子の世界の相転移現象〜」（北沢正清、国広悌二, ISBN:978-4320035492）  
の挿入図生成に用いた数値計算コードの説明と、これらのコードおよび数値計算全般に関するヒントの走り書きである。

コードは  
https://github.com/MasakiyoK/Saizensen/  
に置いてある。

<br />

# コード・実行環境

これらのコードはPython3を使って書かれている。コードの実行には、Python3およびパッケージnumpy, scipy, matplotlib (図2.7のみ、更にcartpy)が必要である。本書執筆時点において、WindowsおよびMacOS上で実行環境を整備するにはanacondaをインストールするのが容易と思われる。anacondaを使うと、Pythonと共にこれらのパッケージもインストールされる。Linuxなら、OSをインストールした段階でPython3が入っていることが多いので、pip環境でパッケージ整備するのが良かろう。いずれもanaconda, pipなどのキーワードで検索を行えばインストール方法はすぐに見つかる。Google Colaboratoryなどのオンライン上の実行環境を使うのも手軽でよいかもしれない。

本コード群の各ファイルは、ファイル名冒頭の数字が書籍内の図番号に対応している。各ファイルをPythonで実行すると数値計算が行われ、計算結果が本書の挿入図とほぼ同じpdfファイルとして出力される。

<br />

# Python初心者向けの一般的なヒント

Pythonは「インタプリタ」であり、c言語などのようにコードをコンパイルすることなくそのまま実行する。
このため手軽に実行できる反面、処理速度はコンパイル言語に比べてはるかに劣る。
筆者の感触としては、標準的な処理においてc/c++言語と比べて数十倍程度遅い。

---

Pythonでは、コードのブロック（cなどの言語で{ }で囲まれる部分）を字下げ（インデント）によって表現する。逆に、字下げにはブロックの指定という明確な意味があるので、適当に字下げしてはいけない。

---

Pythonの変数は、基本的にintやdoubleのような「型」を持たない。変数宣言する必要もないし、違う型を代入しても怒られない。

---

ループ処理にはfor命令を使うが、c言語のように変数の値を変化させる方式ではなく、配列（リスト）の要素に対して行われることに注意。
たとえば
`
for i in range(5):
`
としたときの`range(5)`は`(0,1,2,3,4)`という配列と解釈できるもので、変数`i`にこの配列の要素を代入しながらループが回る。詳しくは、「Python　ループ」などのキーワードで検索。

---

「内包表記」と「ラムダ式」は有用な機能で、本コード群でも多用している。以下にも簡単な説明を載せるが、本コード群を理解するうえでも必須の概念なので、意味が分からない読者は理解しておいてほしい。

---

numpyの「スライス」も、特に多次元配列を扱う際に重宝する強力な機能なので要検索。

---

関数の引数を必要なときだけ変更することのできる、「デフォルト引数」も便利な機能である。c++などと違い、Pythonでは変数名を指定してデフォルト引数に値を代入することができ、この場合変数の順番を気にしなくて良いのが便利である。（c++もこの仕様を採用すればよいのにといつも思う。）

---

Pythonを使いこなす秘訣は、分からないことはネット上で検索しまくることである。
数値計算で現れる処理の多くはありふれたもので、過去に誰かが同じ処理を行うコードを書いたことのあるものである。Pythonではそのような処理がパッケージ内の関数として提供されている。
目的の処理を実行する関数がパッケージの中に存在するならばそれを呼び出せば処理が完了するのだから、問題解決への最短ルートはその関数の名前を知ることに始まる。そして、情報はネット上に溢れている。現代のプログラミングでは、アルゴリズムを書き下す能力より、検索エンジンを使いこなす能力のほうが重要である、といっても過言ではないかもしれない。

なお、筆者はプログラミング言語Pythonに関する書籍を読んだことがない。この言語のことは、全てネット上の情報と人づてで学んだ。（このため、本コード群やこのreadmeには誤った内容が含まれているかもしれない。ネット上で得られる情報は常に怪しいものである。）

---

Pythonには、「ループ(for)を書いたら負け」という格言がある。
これは、Pythonではfor文の処理が重いのと、一見for文が必要になる多くの処理がfor文なしでコーディングでき、多くの場合その方が効率性や可読性に優れるためである。特に配列（リスト）に対する多くの処理はnumpyの活用などによってループなしで圧倒的に短く書ける。

もちろんループが必要になることはあるし、ループを使えば簡易な実装ができる局面で痩せ我慢して複雑なことをしても仕方ないのだが、特に初心者は「ループを書いたら負け」を意識しながらコードを書くことでPythonらしいコードが早く書けるようになると思う。

---

numpy, scipy, matplotlibなどのパッケージは、ネット上にマニュアルが公開されている。基本的に英語なので日本人にとっては言語の壁があるものの、パッケージ内の関数や引数の意味、使い方を正しく理解するためにはこれらのマニュアルを参照するのが良い。ソースコードもウェブ上で閲覧できる。特にある程度以上の細かい処理が必要になった場合には、これらの資料を直接参照するのが早い。

<br />

# 本コードを理解するためのヒントと注意

内包表記はリストを生成する便利な機能であり、本コード群では特にmatplotlibでプロットを行う際に多用している。
例として、`func(x)`という関数があったとして、この関数を区間[0,10]でグラフ表示することを考えよう。
このとき、
```
grid = np.linspace( 0. , 10. , 100 ) #区間[0,10]を100等分したnumpyのリストを作る
result = []
for x in grid:
    result.append( func(x) )
plt.plot( grid , result )
```
のように、y座標を表す配列resultをfor文で作ってからplotしても良いのだが、これと同じことが内包表記を使うと
```
grid = np.linspace( 0. , 10. , 100 )
result = [ func(x) for x in grid ] #内包表記
plt.plot( grid , result )
```
で実現できる。さらに、変数resultの定義を省いて
```
grid = np.linspace( 0. , 10. , 100 )
plt.plot( grid , [ func(x) for x in grid ] )
```
でも良い。こうすることで、もとのコードより3行短くなって得した気分が味わえる。本コード群では、基本的にこのスタイルを採用している。

なお、関数func(x)がnumpyの関数(np.sinとか)と基本的な演算（四則演算等）のみで構成されているならば、
```
grid = np.linspace( 0. , 10. , 100 )
plt.plot( grid , func(grid) )
```
で同じことが実現でき、これが実行速度的にもコード短縮の観点からもベストである。

---

簡単な関数を定義する場合には、ラムダ式を使うと便利なことが多い。
例として、
```
func1 = lambda x : x*x
func2 = lambda x , y : x+y
```
と、
```
def func1(x):
    return x*x
def func2( x , y ):
    return x+y
```
は同じ関数func1, func2を生成する。
（ただし、ラムダ式は本来、無名関数（名前を持たない一時的な関数）を定義するための機能であり、この例のようにラムダ式を変数に代入するのは非推奨らしい。筆者は気にせずに使っているが。）

---

本コード群では積分(scipy.integrate)や最小値探索(scipy.optimize.minimize_scalar)などの、scipyの関数を多用している。
これらの関数の戻り値は積分結果などの解となる値そのものではなく、解の値を含むタプル（変数のリスト）や構造体である。
例えばscipy.integrateの戻り値は積分値、推定誤差などからなるタプルなので、
```
answer = scipy.integrate.quad( func , min , max )
```
としたときにanswerに代入されるのは積分結果の値ではなく、タプルである。answerとして積分値がほしければ、タプルの中の第一変数が積分値なので、
```
answer = scipy.integrate.quad( func , min , max )
answer = answer[0]
```
あるいはさらに省略して、
```
answer = scipy.integrate.quad( func , min , max )[0]
```
とする必要がある。
戻り値がタプルなのは一見複雑な仕様だが、これによって必要に応じて数値積分の推定誤差などの情報にも容易にアクセスできるようになっているのである。

minimize_scalarのときは同様の理由で、
```
answer = scipy.optimize.minimize_scalar( func , ( min , max ) ).x
```
とすることで、answerに極小値の引数xの値が代入される。

---

積分や最小値探索などでエラーメッセージや警告が出る場合には、きちんとメッセージを読んで原因を究明しよう。メッセージの中に、思わぬバグや高速化のヒントが隠されていることが多いためである。特に、積分の収束性が悪いという警告が出たときは、そもそも積分不可能だったり、積分区間に発散があるなどの理由で数値積分が難しい関数である可能性がある。後者の場合には、積分区間を複数に区切る等の工夫によって劇的な高速化・安定化が実現できることがよくある。筆者は、大学院生からポスドク時代に数値積分の高速化が必要な問題に随分悩まされた関係で、数値積分は数値計算の中でも技術力と物理への理解の差が出る腕の見せ所だと思っている。

ただし本コード群ではコードの見やすさを優先する目的でそのような工夫は一切行っていない。特に図3-1と図3-2のHRG熱力学の計算で行う数値積分はそのような工夫を行うべき箇所で、実際これらのコードでは警告が出るしかなりの計算時間を要するのだが、涙を飲んでコードのシンプルさを優先した。また、本コード群では積分ルーチンとしてもっぱら`scipy.integrate.quad`のみを使っているが、他のルーチンを使うことで高速化が見込めることもあるので、研究現場ではそのような検討をしてほしい。

---

本コード群に含まれているファイル plt_setting.py は、matplotlibのプロット関連で筆者が普段使っている設定をまとめたもので、コード冒頭で呼び出すことでこれらの設定がプロット出力に反映される。筆者が原著論文の挿入図を作成するために作り、2022年時点で実際に研究現場で使っているものなので、読者の皆さんが学術利用のための図を作る際にも参考になると思う。

なお、このコード内で行っている配色の変更（35-43行目）については
https://maskit.hatenablog.com/entry/2022/05/29/124006  
が参考になると思う。

x,y軸の目盛間隔を変更する関数xtics, yticsもとても便利。（関数名はgnuplotから拝借。）

---

本コード群では、基本的に各コードで数値計算とグラフ出力をセットで行っているが、これはコード公開の都合上そうしているためで、このようなスタイルが数値計算において常にお勧めできるわけではない。
例として、グラフ描画の微調整を行うような場合には、都度同じ数値計算を繰り返すと処理時間が掛かってしまうので、計算結果を一旦ファイルに出力した上で別のコードでグラフ出力を調整するのがよい。

---

浮動小数点変数には、必ず小数点を付けて整数と区別するクセを付けることをおすすめしたい。実はPythonだと多分問題ないのだが、c/c++など他の言語ではこの区別がバグ防止に役立つ。なお、浮動小数点変数`5.0`はゼロを省いて`5.`と省略でき、筆者はこの記法を好んでいるが、多分一般的ではない。

<br />

# その他数値計算全般に関するヒント

本コード群ではもっぱらPythonを使って数値計算を行っているが、当然ながら全ての数値計算でPythonを使うことが勧められるわけではない。
特に、上でも書いたとおりPythonには処理速度がコンパイラ系の言語と比べて圧倒的に遅いというデメリットがある。このため、見ている間に計算が終了する程度の処理ならともかく、計算時間が数時間を超えるような場合にはc/c++などの処理速度の速い言語の使用を検討したい。その一方で、Pythonにはコード開発が容易で開発時間が短縮できるという大きな利点があり、実行時間がそれほど長くない場合にはこの長所がいかされる。

数値計算の最終目標は科学的な結果を得ることなので、この目的に最短時間で到達できるよう、目的に適した手段を選択することが重要であり、そのためのセンスを磨く必要がある。Mathematicaなども含めた多様な数値計算の手段を選択肢として身に付けておくことが望ましいと思う。

---

上にも書いたとおりPythonの強力な利点は豊富なパッケージであり、Pythonプログラミングではこれらパッケージを最大限に活用することになる。上ではそのようなスタイルを勧めたが、パッケージをブラックボックスとして使っていると思わぬ勘違いから間違った結果を得てしまうことになりがちなので、この点には注意したい。

科学数値計算で最も犯してはならないミスは、間違った計算結果を世間に発表することである。現代のプログラミングにおいて既存のパーツの活用が必須なのは確かだが、この致命的なミスを避けるためには、パッケージを完全なブラックボックスとして扱うのではなく、中身を理解することも常に意識しておきたい。パッケージ内部の理解は、問題解決に適したアルゴリズムを選択し、処理を高速かつ安全に遂行する上でも著しく有用である。

なお、本コード群で行っている積分、最小値探索、非線形方程式の数値解法などのアルゴリズムを解説した日本語文献としては、少し古いがNumerical Recipesをおすすめしたい。

---

思わぬバグを避けるためには、挙動確認のためのテスト計算の手間を惜しまないことも重要である。何かの処理を書いたら、パラメータ依存性などを調べて自然な結果が得られているかを調べるクセをつけよう。特に、物理的な関数ではパラメータの極限的振る舞いの確認が有効である。このような確認を行うことで、課題そのものに対する理解が深まることもよくある。

---

最後に、筆者の経験では、物理の数値計算においてバグやミスを避けるための最も重要かつ有効な秘訣は、計算結果が物理的に妥当かを徹底的に吟味・考察することである。コード上のデバッグも当然重要だが、最後にものをいうのは物理的直感であることを強調してこのreadmeを終えたい。

