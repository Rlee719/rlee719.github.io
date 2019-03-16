---
title: Matrix Analysis
mathjax: true
---

A brief review of Matrix Analysis by R. A. Horn.

<!--more-->

# Chapter 1: Prerequisite
The purpose of this chapter is to introduce some helpful concept and proposition without demonstration, many of which establish foundation for the rest chapter of this book.  

## 0.1 Vector Space

### _0.1.1-0.1.2 Vector Space_

#### _Definition_

**Field**: A set of scalar, denoted as $f$.    

**Vector**: An n-element array formed by element from corresponding filed $f$. Denoted as $V$.  

**Vector Space**: A set of vectors which are closed under addition and multiplication operation, and in addition contain a zero vector.  

### _0.1.3 Subspace & Span_  

#### _Definition_

**Nontrivial Subspace**: Subspace except $\\{0\\}$, $V$.  

**Proper Subspace**: Subspace except $V$.  
   
#### _Generating Subspace_

**Definition**: The intersection of all subspace of $V$ that contain $S$ is the generating subspace of $V$, denoted as $spanS$. The generating subspace of null vector set $span\\{\\}={0}$.  

**Proposition**: $spanS$ equals to all linear combination of vectors in $S$.  

**Sum of Space** : The sum of space $S_1$ and $S_2$:  

$$ S_1 + S_2 = span\{S_1\cup S_2\} $$  

### _0.1.5-0.1.6 Basis_

#### _Definition_

A list of linearly independent vectors from vector space $V$, which take $V$ as the generating subspace of them is a basis of of $V$.  

#### _Proposition_ 

* A list of vectors which generate $V$ are one basis of $V$, only when any of its proper subset can no generate $V$.  
* A list of linearly independent vectors are one basis of $V$ only when any of the proper subset of $V$ that contains the list can not be linearly independent.  
* Every element in $V$ can be represented by basis in a **sole** way.  
* Any list of linearly independent vector in $V$ can be expand to a basis of &V& in some way. (Expand refer to "add one or more vectors to the list".)  

### _0.1.7 Dimension_

#### _Definition_

The amount of vectors in every basis of $V$ is the dimension of $V$, denoting as $dim V$.  
For infinite-dimensional space $V$, there exist a  one-to-one relationship for elements in any two basis of $V$. 
.
#### _Proposition_

$$
\begin{split}
dim \mathbf{R}^n &= n\\
\\
dim \mathbf{C}^n &= \begin{cases} n &\mathbf{F}=\mathbf{S}\\
2n &\mathbf{F}=\mathbf{R}\end{cases}
\end{split}
$$

$R$, $C$ and $S$ refer respectively to real number field, complex number field and pure complex number field.  

#### _Subspace intersection lemma_

$$
dim(S_1 \cap S_2) + dim(S_1 + S_2) = dimS_1 + dimS_2
\nonumber
$$

$Then$  

$$
dim(S_1 \cup S_2) = dimS_1 + simS_2 - dim(S_1 + S_2) \geq dimS_1 + dim S_2 - dimV
$$

### _0.1.8 Isomorphism_

#### _Definition_

$U$, $V$ are vector space in field $F$,  
$F: U \to V$ is an invertible function, that for $x, y \in V, a,b \in F:$

$$
f(ax+by) = af(x) + bf(y)
$$ 
 
Then $f$ is a isomorphism; $U$, $V$ are isomorphic.

#### _Proposition_ 

Two vector space in one field are isomorphic only when they have same dimension

#### _Example of building a $f$_

$If$  
$V$ is a n-dimension space in $F$, $\beta={x_1,\cdots ,x_n}$, $x=a_1 x_1+ \cdots +a_n x_n$, $a_i \in F$, $x_i \in V$  

$Define$  

$$\lbrack x\rbrack _\beta = \lbrack a_1, \cdots , a_n \rbrack ^T$$  

$Then$   
for arbitrary $\beta$, $f:x \to \lbrack x \rbrack _\beta$ is one isomorphic.

## 0.2 Matrix

Denote $A=\lbrack a_{ij} \rbrack \in M_{m,n} (F) $ as a $m \times n$ matrix in field $F$; while $a_{ij}$ denote its element.  

### _0.2.2-0.2.3 Linear transformation_ 

#### _Definition_

$If$  
$U$, $V$ respectively are n-d space and m-d space in $F$, $\beta _U$, $\beta _V$ are basis of $U$, $V$;  
$T$ is a linear transformation that maps $x$, $y$ in $U$, $V$ to $\lbrack x \rbrack _{\beta _U}$ as n-d vectors, $\lbrack y \rbrack _{\beta _V}$ as m-d vectors in $F$.  

$Then$  
construct a $T(x)=y$, which lead to $\lbrack y \rbrack _{\beta V}=A \lbrack x \rbrack _{\beta _U}$.

$A$ is denote as the representation matrix of $T$.  
$T$(or corresponding $A$) depend on $\beta _U$, $\beta _V$.  


**Domain**: $F_n$  

**range**: $rangeA=\\{y \in F^m : y=\mathbf{A}x \\}$, $x \in F^n$  

**null space & Nullity**: $nullspaceA=\\{x \in F^n: \mathbf{A}x=0\\},\ nullityA = dim\ nullspaceA$  

**rank**: $rankA = dim\ range A$

#### _Rank-Nullity theorem_
$$
dim\ rangeA + dim\ nullspaceA = rankA + nullityA =n
$$

### _0.2.5 Conjugate transpose & trace_

#### _Definition_ 

$A^{*} = (\bar{A})^T,\ A \in M_{m,n}(\mathbf{C})$ is the conjugate transpose(Hermitian adjoint) of $A$.  

<br/>

If $A^T=-A$, then $A$ is skew symmetric.  

If $A^{*} A=I $, then $A$ is unitary.  

#### _Division_

for any $A \in M_n (\mathbf{F})$, ${\exists}\ S(A),\ C(A)$ that  
$S(A)=\frac{1}{2} (A + A^T)$ is symmetric.  
$C(A)=\frac{1}{2} (A - A^T)$ is skew symmetric.  
$A = S(A) + C(A)$

#### _Trace_ 

$ tr A = \sum_{k=1}^q a_{kk},\ A=\lbrack a_{ij} \rbrack \in M_{m,n} (\mathbf{F}),\ q=min\\{m,n\\}$  

### 0.2.7 Column Space & Row Space

$rangeA=spanS$ is denoted as the column space of $A$, while $S$ is the list of column vectors of $A$.  

## 0.4 Rank

### _0.4.1 Definition_

For $A \in  M_{m,n}(\mathbf{F})$, $rank A =dim\ range A$ represent the length of the longest linearly independent column or row vector list.  

$$rank A^T = A$$  

### _0.4.2 Linear equations_

#### _Solvability_

If equation $A x = b$ is solvable, then it is consistent.  

$A x = b$ is consistent only when:  

$$
rank \lbrack {A}\quad b \rbrack = rank A
$$
while $\lbrack A\quad b \rbrack$ refers to augmented matrix.  

If $rank\lbrack A\quad b \rbrack = rank A$, then  

$$
dim\ rangeA = dim range \lbrack A\quad b \rbrack\\
dim\ span S_{A_{column}} = dim\ spanS _{\lbrack A\quad b \rbrack _{column}}  
$$

then  

$$
b = \sum_{k=1}^n x_k v_k
$$

$v_k$ denote column vector of $A$ and $\lbrack x_1,\dots,x_k \rbrack = x$.

#### _Proposition_

If $rank A=k,\ A \in M_{m,n}(\mathbf{F})$, than $k=n-dim\ nullspace A$

## 0.5 Nonsingularity

### 0.5.1 Definition 

A linear transformation $T$ or corresponding representation matrix $A$ is nonsingular when it only output $0$ when the input equals to 0.

### 0.5.2 Proposition

* $A \in M_{m,n}(\mathbf{F})$ is singular when $m<n$.  
* $A \in M_{m,n}(\mathbf{F})$ is nonsingular when:  
    * $A^{-1}$ exist.
    * $rankA=n$.
    * $detA \ne 0$.

---
