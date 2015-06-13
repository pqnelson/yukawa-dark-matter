numeric u;
u:=0.5pc;
outputtemplate := "%j-%c.mps";
beginfig(0)
  % axes
  drawdblarrow (-u,0)--(25u,0);
  drawdblarrow (0,-10u)--(0,10u);

  draw (20u,-10u)--(20u,10u) dashed evenly;
  

  draw (0,-7u)--(20u,7u);

  label.rt(btex Initial Guess $m(0)$ etex, (25u,0));
  label.bot(btex Residual etex, (0,-10u));
  dotlabel.lrt(btex $m_{\chi}$ etex, (20u, 0u));
endfig;
beginfig(1)
  % axes
  drawdblarrow (-u,0)--(25u,0);
  drawdblarrow (0,-4u)--(0,4u);

  draw (20u,-4u)--(20u,4u) dashed evenly;
  
  draw ((0,2u)
      for i=0 upto 20u:
      --(i,2u*cosd(i*10))
    endfor );

  label.rt(btex Initial Guess $m(0)$ etex, (25u,0));
  label.bot(btex Residual etex, (0,-4u));
  dotlabel.lrt(btex $m_{\chi}$ etex, (20u, 0u));
endfig;
end;
