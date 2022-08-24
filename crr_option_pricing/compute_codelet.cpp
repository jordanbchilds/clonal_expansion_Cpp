// Copyright (c) 2022 Graphcore Ltd. All rights reserved.
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//      http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// OPTION PRICING: A SIMPLIFIED APPROACH - John Cox, Stephen Ross, Mark Rubinstein - 1979

#include <poplar/Vertex.hpp>
#include <cmath>
#include <algorithm>
#include <array>
#include <ipu_vector_math>

using namespace poplar;

class CRRVertex : public poplar::Vertex
{
public:
  // Fields
  poplar::Input<int> n;
  poplar::Input<float> S;
  poplar::Input<float> K;
  poplar::Input<float> r;
  poplar::Input<float> q;
  poplar::Input<float> v;
  poplar::Input<float> T;
  poplar::Input<bool> isCallOption;
  poplar::Input<bool> isAmerican;

  poplar::InOut<poplar::Vector<float>> prices;

  poplar::Output<float> out;

  float spotCalc(int i, int j, float S, float u, float d)
  {
    // IPU has instructions f32exp and f32log, but there is no f32pow, so the second line is much faster
    // return S * ipu::pow(u, j - i) * ipu::pow(d, i);
    return S * ipu::exp(float(j - i) * ipu::log(u)) * ipu::exp(float(i) * ipu::log(d));
  };

  bool compute()
  {
    int i, j;
    int arraySize = n + 1;

    float dt = T / n;
    float u = ipu::exp(v * ipu::sqrt(dt));
    float d = 1 / u;
    float p = (ipu::exp((r - q) * dt) - d) / (u - d);
    float exprdt = ipu::exp(-r * dt);

    // build the array of option valuations at expiration date
    if (isCallOption)
    {
      for (i = 0; i <= n; i++)
      {
        prices[i] = ipu::fmax(spotCalc(i, n, S, u, d) - K, float(0.0));
      }
    }
    else
    {
      for (i = 0; i <= n; i++)
      {
        prices[i] = ipu::fmax(K - spotCalc(i, n, S, u, d), float(0.0));
      }
    }

    // iterate back through by calculating binomial value
    if (!isAmerican)
    {
      for (j = n - 1; j >= 0; j--)
      {
        for (i = 0; i <= j; i++)
        {
          prices[i] = exprdt * (p * (prices[i]) +
                                           (1 - p) * (prices[i + 1]));
        }
      }
    }
    else if (isCallOption)
    {
      for (j = n - 1; j >= 0; j--)
      {
        for (i = 0; i <= j; i++)
        {
          prices[i] = ipu::fmax(spotCalc(i, j, S, u, d) - K,
                               exprdt *
                                   (p * (prices[i]) +
                                    (1 - p) * (prices[i + 1])));
        }
      }
    }
    else
    {
      for (j = n - 1; j >= 0; j--)
      {
        for (i = 0; i <= j; i++)
        {
          prices[i] = ipu::fmax(K - spotCalc(i, j, S, u, d),
                               exprdt *
                                   (p * (prices[i]) +
                                    (1 - p) * (prices[i + 1])));
        }
      }
    }

    *out = prices[0];
    return true;
  }
};
