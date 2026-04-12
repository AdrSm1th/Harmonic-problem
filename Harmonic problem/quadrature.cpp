//quadrature.cpp

#include "quadrature.h"

std::vector<double> Quadrature::getPoints(int order) {
   switch (order) {
   case 1: return { 0.0 };
   case 2: return { -0.5773502691896257, 0.5773502691896257 };
   case 3: return { -0.7745966692414834, 0.0, 0.7745966692414834 };
   }
}

std::vector<double> Quadrature::getWeights(int order) {
   switch (order) {
   case 1: return { 2.0 };
   case 2: return { 1.0, 1.0 };
   case 3: return { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556 };
   }
}