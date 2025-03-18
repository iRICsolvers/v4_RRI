# RRI on iRIC

# Japanese

## RRI on iRICについて

RRI on iRICは、土木研究所ICHARMが開発したRRIモデルをiRICで使えるよう、RIverlinkの旭氏が中心となって、RRIモデルを改変したプログラムです。その後、土木研究所ICHARMが開発している降雨-土砂流出モデル等を統合したものです。
iRICのポストプロセッサを使って計算結果の表示を行うなど、便利に利用できるようになりました。

## 動作環境

iRICに同封されているものは，Intel Open APIを使ってコンパイルをしています．

同環境で使用したい場合は例えば下記のページなどで環境構築をお願いします．

https://qiita.com/Kazutake/items/a069f86d21ca43b6c153

## 改良版をiRICで動かす

make.batでコンパイルできます．もちろん，そのほかの方法でコンパイルいただいても大丈夫だと思います．

コンパイルするとrri.exeが生成されるので，iRICのインストールフォルダにある，solvers/RRIにコピーすると更新ソルバーが使えます．

条件設定ファイルdefinition.xmlを変更した際も同様です．

# English

## About RRI on iRIC

RRI on iRIC is a modified version of the RRI model, originally developed by ICHARM (International Centre for Water Hazard and Risk Management) at PWRI (Public Works Research Institute). This adaptation was primarily led by Mr. Asahi of RIverlink to make the RRI model compatible with iRIC. Subsequently, it has been integrated with rainfall-sediment runoff models and other functionalities developed by ICHARM. It now offers enhanced usability, including the display of calculation results using iRIC's post-processor.

## Operating Environment

The version packaged with iRIC is compiled using the Intel Open API.
If you wish to use it in the same environment, please follow the instructions on a resource like the following page for environment setup:
https://qiita.com/Kazutake/items/a069f86d21ca43b6c153

## Running the Improved Version with iRIC

You can compile it using make.bat. Of course, you may also compile it using other methods.
Compiling generates rri.exe. Copying this to the solvers/RRI directory within your iRIC installation folder will update the solver.
The same procedure applies when you modify the definition file definition.xml.

