% Basic setup
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{evanbook}

\LoadClass[a4paper, twocolumn]{book}

%\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}

%[style=verbose]
\usepackage[columnsep=.5cm, top=2cm, bottom=2cm, left=1cm, right=1cm]{geometry}
\usepackage{makeidx}
\usepackage{import}
\usepackage{bigints}
\usepackage{algorithm}
\usepackage{algpseudocode}
\usepackage{mathtools}
\usepackage{amsmath,amsthm,amssymb,amscd}
\usepackage{graphicx}
\usepackage{float}
\usepackage{multicol}
\usepackage{tabularx}
\usepackage[format=plain,
            labelfont={bf},
            textfont=it]{caption}
\usepackage{subcaption}
\usepackage[nice]{units}
\usepackage[mathscr]{euscript}
\usepackage{stmaryrd}
%\usepackage{fouriernc}
%\usepackage{igsupmath}
\usepackage{amsbsy}
\usepackage{amssymb}
\usepackage{bm}
\usepackage{xspace}
\usepackage[scaled=0.90]{beramono}
\usepackage[10pt]{moresize}
\usepackage{makecell}
%\usepackage[nameinlink, capitalize, noabbrev]{cleveref}
%\crefformat{equation}{(#2#1#3)}
%\crefrangeformat{equation}{(#3#1#4)--(#5#2#6)}
%\crefmultiformat{equation}{(#2#1#3)}%
%{ and~(#2#1#3)}{, (#2#1#3)}{ and~(#2#1#3)}

% url links :
\usepackage{hyperref}

\hypersetup{
  final,
  colorlinks,
  linkcolor={red!80!black},
  citecolor={red!80!black},
  urlcolor={blue!60!red}
}

% right cases :
%\newenvironment{rcases}
%  {\left.\begin{aligned}}
%  {\end{aligned}\right\rbrace}

% fonts :
\DeclareMathAlphabet{\mathcal}{OMS}{cmsy}{m}{n}   % or txsy
\SetMathAlphabet{\mathcal}{bold}{OMS}{cmsy}{b}{n} % or txsy
\DeclareMathAlphabet{\pazocal}{OMS}{zplm}{m}{n}

\newtheoremstyle{remark}%      % Name
  {}%                          % Space above
  {}%                          % Space below
  {}%                          % Body font
  {}%                          % Indent amount
  {\bfseries}%                 % Theorem head font
  {.}%                         % Punctuation after theorem head
  { }%                         % Space after theorem head, ' ', or \newline
  {}%                          % Theorem head spec (can be left empty, meaning `normal')

\newtheoremstyle{definition}%  % Name
  {}%                          % Space above
  {}%                          % Space below
  {}%                          % Body font
  {}%                          % Indent amount
  {\bfseries}%                 % Theorem head font
  {.}%                         % Punctuation after theorem head
  { }%                         % Space after theorem head, ' ', or \newline
  {}%                          % Theorem head spec (can be left empty, meaning `normal')

\newtheorem{theorem}{Theorem}[section]
\newtheorem{corollary}{Corollary}[theorem]
\newtheorem{lemma}[theorem]{Lemma}
\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}

\newcommand{\N}{\mathbb{N}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\bigo}{\mathcal{O}}
\newcommand{\G}{\mathcal{G}}
\newcommand{\V}{\mathcal{V}}
\newcommand{\E}{\mathcal{E}}
\newcommand{\K}{\mathcal{K}}
\newcommand{\Hu}{\mathcal{H}}
\newcommand{\Lu}{\mathcal{L}}
\newcommand{\T}{^\intercal}
\newcommand{\var}{\mathrm{Var}}
\newcommand{\mean}{\mathrm{mean}}
\newcommand{\Ex}{\mathrm{E}}
\newcommand{\tr}{\mathrm{tr}}
\newcommand{\diag}{\mathrm{diag}}
\newcommand{\ones}{\mathrm{ones}}
\newcommand{\mat}{\mathbf}
\newcommand{\CSLVR}{\href{http://cslvr.readthedocs.io/en/latest/}{CSLVR}\xspace}

% some commands for sub- and super-scripts :
\newcommand{\fint}{\rankone{f}^{\mathrm{int}}}
\newcommand{\fext}{\rankone{f}^{\mathrm{ext}}}
\newcommand{\nthd}{N_{\mathrm{thd}}}
\newcommand{\p}{\mathrm{p}}
\newcommand{\pbas}{\varphi}
\newcommand{\ebas}{\phi}
\newcommand{\extrap}{\widetilde}
\newcommand{\melt}{\mathrm{m}}
\newcommand{\normal}{\rankone{\hat{n}}}

%===============================================================================
% variable subscripts :
\newcommand{\ice}{\mathrm{i}} 
\newcommand{\water}{\mathrm{w}}
\newcommand{\liquid}{\mathrm{\ell}}
\newcommand{\solid}{\mathrm{s}}
\newcommand{\vapor}{\mathrm{v}}
\newcommand{\rock}{\mathrm{r}}
\newcommand{\seawater}{\mathrm{sw}}
\newcommand{\sensible}{\mathrm{s}}
\newcommand{\latent}{\mathrm{l}}
\newcommand{\ini}{\mathrm{0}}
\newcommand{\crit}{\mathrm{c}}
\newcommand{\m}{\mathrm{m}}
\newcommand{\e}{\mathrm{e}}
\newcommand{\air}{\mathrm{A}}
\newcommand{\wat}{\mathrm{W}}
\newcommand{\bed}{\mathrm{B}}
\newcommand{\srf}{\mathrm{S}}
\newcommand{\C}{\mathrm{C}}
\newcommand{\tmp}{\mathrm{T}}
\newcommand{\gnd}{\mathrm{G}}
\newcommand{\ext}{\mathrm{E}}
\newcommand{\lat}{\mathrm{L}}
\newcommand{\FS}{\mathrm{FS}}
\newcommand{\BP}{\mathrm{BP}}
\newcommand{\PS}{\mathrm{PS}}
\newcommand{\RS}{\mathrm{RS}}
\newcommand{\IE}{\mathrm{IE}}
\newcommand{\BV}{\mathrm{BV}}
\newcommand{\AG}{\mathrm{AG}}
\newcommand{\geo}{\mathrm{geo}}
\newcommand{\fric}{\mathrm{fric}}
\newcommand{\ob}{\mathrm{ob}}
\newcommand{\eq}{\mathrm{eq}}
\newcommand{\extern}{\mathrm{ext}}
\newcommand{\intern}{\mathrm{int}}
\newcommand{\stz}{\mathrm{stz}}
\newcommand{\lin}{\mathrm{lin}}
\newcommand{\body}{\mathrm{b}}
\newcommand{\ana}{\mathrm{a}}
\newcommand{\true}{\mathrm{true}}
\newcommand{\fs}{\mathring{F}} 
\newcommand{\smb}{\mathring{S}} 
\newcommand{\bmb}{\mathring{B}} 

% italic characters :
\let\Gamma\varGamma
\let\Omega\varOmega
\let\Sigma\varSigma

% Typesetting macros
\newcommand{\blankpage}{\newpage\hbox{}\thispagestyle{empty}\newpage}
\newcommand{\hairsp}{\hspace{1pt}}
\newcommand{\ie}{\textit{i.\hairsp{}e.}\xspace}
\newcommand{\eg}{\textit{e.\hairsp{}g.}\xspace}
\newcommand{\tab}{\hspace*{2em}}
\newcommand{\emp}[1]{\texttt{#1}}
\newcommand{\etc}{\emph{etc}}
\newcommand{\viceversa}{\emph{vice versa}}
\newcommand{\cf}{\textit{cf.}\xspace}

% Math notation macros
\renewcommand{\d}[1]{\ \mathrm{d}#1}
\newcommand{\diff}[1]{\mathrm{d}#1}

\newcommand{\Hdiv}{H(\mathrm{div})}
\newcommand{\Hcurl}{H(\mathrm{curl})}
\newcommand{\Rset}{\ensuremath{\mathbb{R}}}
\newcommand{\inner}[2]{\langle #1, #2 \rangle}
\newcommand{\renni}[2]{\langle #2, #1 \rangle}
\newcommand{\totder}[3][]{\frac{\mathrm{d}^{#1} #2}{\mathrm{d} #3^{#1} }}
\newcommand{\parder}[3][]{\frac{\partial^{#1} #2}{\partial #3^{#1} }}
\newcommand{\delder}[3][]{\frac{\delta{^{#1}} #2}{\delta #3{^{#1}} }}
\newcommand{\deltat}{\Delta t}
%\newcommand{\bm}[1]{#1}
\newcommand{\brac}[1]{\left( {#1} \right)}
\newcommand{\bracc}[1]{\left\{ {#1} \right\}}
\newcommand{\foralls}{\forall\,}
\newcommand{\mesh}{\mathcal{T}_h}
\newcommand{\Poly}[1]{\mathcal{P}_{#1}}
\newcommand{\CiarletSpace}{\mathcal{V}}
\newcommand{\avg}[1]{\langle {#1} \rangle}
\newcommand{\jump}[1]{\left\llbracket {#1} \right\rrbracket}
\newcommand{\jumpn}[1]{\llbracket {#1} \rrbracket_n}
\newcommand{\cgone}{$\mathrm{cG}(1)$}
\newcommand{\supp}{\mathrm{supp}}
\newcommand{\coefvec}[1]{\bm{#1}}
\newcommand{\rankone}[1]{\underline{#1}}
\newcommand{\ranktwo}[1]{\underline{\underline{#1}}}
\newcommand{\rankfour}[1]{\underline{\underline{\underline{\underline{#1}}}}}

% Math operatorsg
\DeclareMathOperator{\Div}{div}
\DeclareMathOperator{\Grad}{grad}

% Difficult names
\newcommand{\nedelec}{N\'ed\'elec}
\newcommand{\babuska}{Babu\v{s}ka}

% Misc macros
\newcommand{\apriori}{\emph{a~priori}}
\newcommand{\apost}{\emph{a~posteriori}}

\newcommand{\testspace}{S_0^h \subset \mathcal{H}_{E_0}^1(\Omega)}
\newcommand{\trialspace}{S_E^h \subset \mathcal{H}_E^1(\Omega)}
\newcommand{\ltwospace}{M^h \subset L^2(\Omega)}

\newcommand{\Q}{\mathcal{Q}}
  
\DeclareMathOperator*{\argmin}{argmin}
\DeclareMathOperator*{\argmax}{argmax}
\newcommand*{\argminl}{\argmin\limits}
\newcommand*{\argmaxl}{\argmax\limits}

\newfloat{myalgo}{tbhp}{mya}

\newenvironment{Algorithm}[2][tbh]%
{\begin{myalgo}[#1]
\centering
\begin{minipage}{#2}
\begin{algorithm}[H]}%
{\end{algorithm}
\end{minipage}
\end{myalgo}}

\def\changemargin#1#2{\list{}{\rightmargin#2\leftmargin#1}\item[]}
\let\endchangemargin=\endlist 

\newcommand{\sups}[1]{\ensuremath{^{\textrm{#1}}}}
\newcommand{\subs}[1]{\ensuremath{_{\textrm{#1}}}}

\newcommand{\specialcell}[2][c]
{
  \begin{tabular}[#1]{@{}c@{}}#2\end{tabular}
}

\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

\makeatletter
\newsavebox{\mybox}\newsavebox{\mysim}
\newcommand{\distras}[1]
{
  \savebox{\mybox}{\hbox{\kern3pt$\scriptstyle#1$\kern3pt}}%
  \savebox{\mysim}{\hbox{$\sim$}}%
  \mathbin{\overset{#1}{\kern\z@\resizebox{\wd\mybox}{\ht\mysim}{$\sim$}}}%
}
\makeatother

% aligned matrix negative signs:
\makeatletter
\renewcommand*\env@matrix[1][c]{\hskip -\arraycolsep
  \let\@ifnextchar\new@ifnextchar
  \array{*\c@MaxMatrixCols #1}}
\makeatother

%===============================================================================
% code highlighting :
\usepackage[final]{listings}

% define custom colors :
\usepackage{color}
\usepackage{xcolor}

\definecolor{bg}{rgb}{0.96,0.96,0.85}
\definecolor{deepblue}{rgb}{0,0,0.5}
\definecolor{deepred}{rgb}{0.6,0,0}
\definecolor{deepgreen}{rgb}{0,0.5,0}
\definecolor{gray}{gray}{0.5}

\renewcommand{\lstlistlistingname}{Code Listings}
\renewcommand{\lstlistingname}{Code Listing}


\colorlet{commentcolour}{cyan!40!black}
\colorlet{stringcolour}{magenta!40!blue}
\colorlet{keywordcolour}{red!90!black}
\colorlet{exceptioncolour}{yellow!50!red}
\colorlet{commandcolour}{blue!60!black}
\colorlet{numpycolour}{green!70!blue!90!black}
\colorlet{digitcolour}{cyan!70!blue!90!black}
\colorlet{literatecolour}{magenta!90!black}
\colorlet{promptcolour}{green!80!black}
\colorlet{specmethodcolour}{violet}
\colorlet{indendifiercolour}{green!70!white}

\newcommand{\framemargin}{5ex}

\newcommand{\literatecolour}{\textcolor{literatecolour}}
\newcommand{\digitcolour}{\textcolor{digitcolour}}

\newcommand\Small{\fontsize{1.00}{5.0}\selectfont}

\newcommand\pythonstyle{\lstset{
%keepspaces=true,
language=python,
showtabs=true,
tab=,
tabsize=2,
basicstyle=\ttfamily\ssmall,%\setstretch{.5},
stringstyle=\color{stringcolour},
showstringspaces=false,
alsoletter={1234567890},
otherkeywords={\ , \}, \{, \%, \&, \|},
keywordstyle=\color{keywordcolour}\bfseries,
emph={and,break,class,continue,def,yield,del,elif ,else,%
except,exec,finally,for,from,global,if,import,in,%
lambda,not,or,pass,print,raise,return,try,while,assert},
emphstyle=\color{blue}\bfseries,
emph={[2]True, False, None},
emphstyle=[2]\color{keywordcolour},
emph={[3]object,type,isinstance,copy,deepcopy,zip,enumerate,reversed,list,len,dict,tuple,xrange,append,execfile,real,imag,reduce,str,repr},
emphstyle=[3]\color{commandcolour},
emph={Exception,NameError,IndexError,SyntaxError,TypeError,ValueError,OverflowError,ZeroDivisionError},
emphstyle=\color{exceptioncolour}\bfseries,
%upquote=true,
morestring=[s]{"""}{"""},
morestring=[s]{'''}{'''},
commentstyle=\color{commentcolour}\slshape,
%emph={[4]1, 2, 3, 4, 5, 6, 7, 8, 9, 0},
emph={[4]ode, fsolve, sqrt, exp, sin, cos, arccos, pi,  array, norm, solve, dot, arange, , isscalar, max, sum, flatten, shape, reshape, find, any, all, abs, linspace, legend, quad, polyval,polyfit, hstack, concatenate,vstack,column_stack,empty,zeros,ones,rand,vander,grid,pcolor,eig,eigs,eigvals,svd,qr,tan,det,logspace,roll,min,mean,cumsum,cumprod,diff,vectorize,lstsq,cla,eye,xlabel,ylabel,squeeze,plot,median,std,hist},
emphstyle=[4]\color{numpycolour},
emph={[5]__init__,__add__,__mul__,__div__,__sub__,__call__,__getitem__,__setitem__,__eq__,__ne__,__nonzero__,__rmul__,__radd__,__repr__,__str__,__get__,__truediv__,__pow__,__name__,__future__,__all__},
emphstyle=[5]\color{specmethodcolour},
emph={[6]assert,range,yield},
emphstyle=[6]\color{keywordcolour}\bfseries,
% emph={[7]self},
% emphstyle=[7]\bfseries,
literate=*%
{:}{{\literatecolour:}}{1}%
{=}{{\literatecolour=}}{1}%
{-}{{\literatecolour-}}{1}%
{+}{{\literatecolour+}}{1}%
{*}{{\literatecolour*}}{1}%
{/}{{\literatecolour/}}{1}%
{!}{{\literatecolour!}}{1}%
%{(}{{\literatecolour(}}{1}%
%{)}{{\literatecolour)}}{1}%
{[}{{\literatecolour[}}{1}%
{]}{{\literatecolour]}}{1}%
{<}{{\literatecolour<}}{1}%
{>}{{\literatecolour>}}{1}%
{>>>}{{\textcolor{promptcolour}{>>>}}}{1}%
%{0}{{\digitcolour{0}}}{1}%
%{1}{{\digitcolour{1}}}{1}%
%{2}{{\digitcolour{2}}}{1}%
%{3}{{\digitcolour{3}}}{1}%
%{4}{{\digitcolour{4}}}{1}%
%{5}{{\digitcolour{5}}}{1}%
%{6}{{\digitcolour{6}}}{1}%
%{7}{{\digitcolour{7}}}{1}%
%{8}{{\digitcolour{8}}}{1}%
%{9}{{\digitcolour{9}}}{1}%
%{.0}{{\digitcolour{.0}}}{2}%
%{.1}{{\digitcolour{.1}}}{2}%
%{.2}{{\digitcolour{.2}}}{2}%
%{.3}{{\digitcolour{.3}}}{2}%
%{.4}{{\digitcolour{.4}}}{2}%
%{.5}{{\digitcolour{.5}}}{2}%
%{.6}{{\digitcolour{.6}}}{2}%
%{.7}{{\digitcolour{.7}}}{2}%
%{.8}{{\digitcolour{.8}}}{2}%
%{.9}{{\digitcolour{.9}}}{2}%
%{e0}{{\digitcolour{e0}}}{2}%
%{e1}{{\digitcolour{e1}}}{2}%
%{e2}{{\digitcolour{e2}}}{2}%
%{e3}{{\digitcolour{e3}}}{2}%
%{e4}{{\digitcolour{e4}}}{2}%
%{e5}{{\digitcolour{e5}}}{2}%
%{e6}{{\digitcolour{e6}}}{2}%
%{e7}{{\digitcolour{e7}}}{2}%
%{e8}{{\digitcolour{e8}}}{2}%
%{e9}{{\digitcolour{e9}}}{2}%
%{e-}{{\digitcolour{e-}}}{2}%
%{\ }{{ }}{1}%
,%
breaklines=true,
breakatwhitespace= true,
%xleftmargin=\framemargin,
%xrightmargin=\framemargin,
aboveskip=1ex,
frame=trbl,
%frameround=tttt,
rulecolor=\color{black!40},
%framexleftmargin=\framemargin,
%framextopmargin=.1ex,
%framexbottommargin=.1ex,
%framexrightmargin=\framemargin,
%framexleftmargin=1mm, framextopmargin=1mm, frame=shadowbox, rulesepcolor=\color{blue},#1
%frame=tb,
backgroundcolor=\color{black!5}
}}

% Python environment
\lstnewenvironment{python}[1][]
{
  \pythonstyle
  \lstset{#1}
}
{}

% Python for external files
\newcommand{\pythonexternal}[2][]
{{
  \pythonstyle
  \lstinputlisting[#1]{#2}
}}

% Python for inline
\newcommand\pythoninline[1]
{{
  \pythonstyle
  \lstinline!#1!
}}

% answer box :
\newcommand\answer[1]
{{
  \begin{center}
    \parbox[t]{5in}{#1}
  \end{center}
}}

% end code highlighting
%===============================================================================


% border matrix with brackets :
\usepackage{etoolbox}
\let\bbordermatrix\bordermatrix
\patchcmd{\bbordermatrix}{8.75}{4.75}{}{}
\patchcmd{\bbordermatrix}{\left(}{\left[}{}{}
\patchcmd{\bbordermatrix}{\right)}{\right]}{}{}

% end function.tex
%===============================================================================

%\usepackage[backend=biber, style=authoryear, sorting=nyt,natbib]{biblatex}
\usepackage[backend      = bibtex,
            maxbibnames  = 2,
            maxcitenames = 2,
            firstinits   = true,
            style        = authoryear,
            sorting      = nyt,
            natbib]{biblatex}
\DefineBibliographyStrings{english}{%
    andothers = {\em et\addabbrvspace al\adddot}
}
\AtEveryBibitem
{
  \clearfield{issn} % Remove issn
  %\clearfield{doi} % Remove doi
  \clearfield{eprint}

  % Remove url except for @online
  \ifentrytype{online}{}
  {
    \clearfield{url}
  }
}

% #1: sorting key, #2: prefix for citation, #3: prefix for bibliography
\DeclareRobustCommand{\VAN}[3]{#2} % set up for citation

\addbibresource{/home/pf4d/cummings.evan@gmail.com/biblio.bib}
