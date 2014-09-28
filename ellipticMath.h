#ifndef ELLIPTICMATH_H
#define ELLIPTICMATH_H
#include <math.h>
#include <gmpxx.h>
struct ec
{
    mpz_class x;//координата x
    mpz_class y;//координата y
};
struct eq
{
    mpz_class a;//первый аргумент уравнения
    mpz_class b;//второй аргумент уравнения
    mpz_class field;//модуль ЭК
    mpz_class q;//порядок циклической подгруппы группы точек ЭК
};

mpz_class invert_element(mpz_class a, mpz_class b)//нахождения обратного элемента
{
    mpz_class field = b;
  mpz_class q, r, x1, x2, y1, y2, d, x, y;
  if (b == 0)
  {
    d = a, x = 1, y = 0;
    return x;
  }
  mpz_class mid;
  x2 = 1, x1 = 0, y2 = 0, y1 = 1;
  while (b > 0) {
    mpz_fdiv_q(q.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    mpz_mul(mid.get_mpz_t(), q.get_mpz_t(), b.get_mpz_t());
    mpz_sub(r.get_mpz_t(), a.get_mpz_t(), mid.get_mpz_t());
    //q = a / b, r = a - q * b;
    mpz_mul(mid.get_mpz_t(), q.get_mpz_t(), x1.get_mpz_t());
    mpz_sub(x.get_mpz_t(), x2.get_mpz_t(), mid.get_mpz_t());
    mpz_mul(mid.get_mpz_t(), q.get_mpz_t(), y1.get_mpz_t());
    mpz_mul(y.get_mpz_t(), y2.get_mpz_t(), mid.get_mpz_t());
    //x = x2 - q * x1, y = y2 - q * y1;

    a = b, b = r;
    x2 = x1, x1 = x, y2 = y1, y1 = y;
  }
  d = a, x = x2, y = y2;

  if (x < 0)
      mpz_add(x.get_mpz_t(), x.get_mpz_t(), field.get_mpz_t());
        //x+= field;
  return x;
};

ec two_ecp(ec p1, ec p2, eq equation)//сложение двух точек
{
      ec p3;
      mpz_class mid = 0;
      mpz_class dy;// = p2.y - p1.y;
      mpz_sub(dy.get_mpz_t(), p2.y.get_mpz_t(), p1.y.get_mpz_t());
      mpz_class dx;// = p2.x - p1.x;
      mpz_sub(dx.get_mpz_t(), p2.x.get_mpz_t(), p1.x.get_mpz_t());

      if (dx < 0)
        //dx+= equation.field;
      mpz_add(dx.get_mpz_t(), dx.get_mpz_t(), equation.field.get_mpz_t());
      if (dy < 0)
        //dy+=equation.field;
      mpz_add(dy.get_mpz_t(), dy.get_mpz_t(), equation.field.get_mpz_t());

    dx=invert_element(dx, equation.field);

      mpz_class l;// = (dy * dx) % equation.field;
      mpz_mul(mid.get_mpz_t(), dy.get_mpz_t(), dx.get_mpz_t());
      mpz_mod(l.get_mpz_t(), mid.get_mpz_t(), equation.field.get_mpz_t());
      if (l < 0)
        //l += equation.field;
      mpz_add(l.get_mpz_t(), l.get_mpz_t(), equation.field.get_mpz_t());
      //p3.x = (l * l - p1.x - p2.x) % equation.field;
      mpz_mul(mid.get_mpz_t(), l.get_mpz_t(), l.get_mpz_t());
      mpz_sub(mid.get_mpz_t(), mid.get_mpz_t(), p1.x.get_mpz_t());
      mpz_sub(mid.get_mpz_t(), mid.get_mpz_t(), p2.x.get_mpz_t());
      mpz_mod(p3.x.get_mpz_t(), mid.get_mpz_t(), equation.field.get_mpz_t());
      //p3.y = (l * (p1.x - p3.x) - p1.y) % equation.field;
      mpz_sub(mid.get_mpz_t(), p1.x.get_mpz_t(), p3.x.get_mpz_t());
      mpz_mul(mid.get_mpz_t(), mid.get_mpz_t(), l.get_mpz_t());
      mpz_sub(mid.get_mpz_t(), mid.get_mpz_t(), p1.y.get_mpz_t());
      mpz_mod(p3.y.get_mpz_t(), mid.get_mpz_t(), equation.field.get_mpz_t());

      if (p3.x < 0)
          mpz_add(p3.x.get_mpz_t(), p3.x.get_mpz_t(), equation.field.get_mpz_t());
        //p3.x += equation.field;
      if (p3.y < 0)
          mpz_add(p3.y.get_mpz_t(), p3.y.get_mpz_t(), equation.field.get_mpz_t());
        //p3.y += equation.field;
     return p3;
};

ec double_ecp(ec p, eq equation)//нахождения удвоенной точки
{
      ec p2;

      mpz_class mid = 0;
      mpz_class dy;// = 3 * p.x * p.x + equation.a;
      mpz_mul_ui(dy.get_mpz_t(), p.x.get_mpz_t(), 3);
      mpz_mul(dy.get_mpz_t(), dy.get_mpz_t(), p.x.get_mpz_t());
      mpz_add(dy.get_mpz_t(), dy.get_mpz_t(), equation.a.get_mpz_t());
      mpz_class dx;// = 2 * p.y;
      mpz_mul_ui(dx.get_mpz_t(), p.y.get_mpz_t(), 2);


      if (dx < 0)
        //dx += equation.field;
      mpz_add(dx.get_mpz_t(), dx.get_mpz_t(), equation.field.get_mpz_t());
      if (dy < 0)
        //dy += equation.field;
      mpz_add(dy.get_mpz_t(), dy.get_mpz_t(), equation.field.get_mpz_t());


       dx=invert_element(dx, equation.field);

      mpz_class l = (dy * dx)%equation.field;
      mpz_mul(mid.get_mpz_t(), dy.get_mpz_t(), dx.get_mpz_t());
      mpz_mod(l.get_mpz_t(), mid.get_mpz_t(), equation.field.get_mpz_t());
      //p2.x = (l*l - p.x - p.x) % equation.field;
      mpz_mul(mid.get_mpz_t(), l.get_mpz_t(), l.get_mpz_t());
      mpz_sub(mid.get_mpz_t(), mid.get_mpz_t(), p.x.get_mpz_t());
      mpz_sub(mid.get_mpz_t(), mid.get_mpz_t(), p.x.get_mpz_t());
      mpz_mod(p2.x.get_mpz_t(), mid.get_mpz_t(), equation.field.get_mpz_t());
      //p2.y = (l*(p.x - p2.x) - p.y) % equation.field;
      mpz_sub(mid.get_mpz_t(), p.x.get_mpz_t(), p2.x.get_mpz_t());
      mpz_mul(mid.get_mpz_t(), mid.get_mpz_t(), l.get_mpz_t());
      mpz_sub(mid.get_mpz_t(), mid.get_mpz_t(), p.y.get_mpz_t());
      mpz_mod(p2.y.get_mpz_t(), mid.get_mpz_t(), equation.field.get_mpz_t());
      if (p2.x < 0)
        //p2.x += equation.field;
      mpz_add(p2.x.get_mpz_t(), p2.x.get_mpz_t(), equation.field.get_mpz_t());
      if (p2.y < 0)
        //p2.y += equation.field;
      mpz_add(p2.y.get_mpz_t(), p2.y.get_mpz_t(), equation.field.get_mpz_t());
      return p2;
};

ec multiply_ecp(ec p, eq equation, mpz_class k)//умноежение точки на число
{
    ec temp;
    temp=p;
    k--;
    while(k!=0)
    {
        if((k%2)!=0)
        {
            if((temp.x==p.x)&&(temp.y==p.y))
                temp=double_ecp(temp, equation);
            else
                temp=two_ecp(temp, p, equation);
            k--;
        }
        k/=2;
        p=double_ecp(p, equation);
    }
    return temp;
};
#endif // ELLIPTICMATH_H
