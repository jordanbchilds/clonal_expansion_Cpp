#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poplar/Graph.hpp>
#include "codelets/random.hpp"

using namespace std;

const size_t workers = 6;  // TODO: should enumerate this
const size_t per_worker = 8;
const size_t random_output_floats = workers * per_worker;

bool use_random = true;
uint32_t random_seed = 0;
uint32_t base = 0, tile_factor = 0, worker_factor = 0;

class random_output_stream_callback : public poplar::StreamCallback {
    public:
        random_output_stream_callback() { }

        void fetch(void *p) override {
            // Just print results.
            float* random_output = static_cast<float*>(p);
            cout << endl << "Results:" << endl;
            for (uint32_t e = 0; e < random_output_floats;) {
              for (uint32_t r = 0; r <  per_worker; ++r) {
#if USE_URAND32
                 cout << static_cast<uint32_t>(random_output[e]) << " ";
#else
                 cout << random_output[e] << " ";
#endif
                 ++e;
              }
              cout << endl;
            }
        }

        poplar::StreamCallback::Result prefetch(void* p) override {
            fetch(p);
            return poplar::StreamCallback::Result::Success;
        }

        void complete() override {};

};

void ipu_program(size_t random_output_floats) {
    const uint32_t num_ipus = 1;
    const uint32_t num_tiles = 1;

    // Enumerate IPUs.
    poplar::DeviceManager manager{};
    auto devices = manager.getDevices(poplar::TargetType::IPU, num_ipus);
    if (devices.empty()) {
        throw poplar::runtime_error("No IPU devices available num_ipus=" + num_ipus);
    }

    auto& target = devices[0].getTarget();

    auto graph = poplar::Graph{target};

    auto random_graph = graph.createVirtualGraph(num_tiles);

    const auto output_element_type = poplar::FLOAT;

    // Create the random output tensor.
    auto random_tensor_out = random_graph.addVariable(
        poplar::FLOAT,
        {random_output_floats},
        "random_tensor_out");

    // Map random I/O tensors to tile(s).
    assert(num_tiles == 1);

    // Program part for setting seeds.
    auto set_seeds1_prog = poplar::program::Sequence{};

    // Add setSeeds op to program,
    if (use_random)
      random_ipu::setRandomSeeds(graph, random_seed, set_seeds1_prog);
    else
      random_ipu::setSeeds(graph, base,tile_factor,worker_factor, set_seeds1_prog);

    // Program part for compute random stream.
    auto random_prog = poplar::program::Sequence{};

    // Add random op to program.
    random_ipu::random(graph, random_tensor_out, random_prog);

    // Output stream decl.
    auto random_output_datastream = random_graph.addDeviceToHostFIFO(
        "random_output_stream",
        poplar::FLOAT,
        random_output_floats);

    // Program parts for I/O.
    auto random_output_prog = poplar::program::Copy(random_tensor_out, random_output_datastream);

    // Process N times.
    auto n_random = poplar::program::Sequence{
        set_seeds1_prog,      // Set seeds1 first
        // First iteration:
        random_prog,          // Random stream
        random_output_prog,   // Output results
#if 0
        // Repeat iterations:
        random_prog,
        random_output_prog,
        random_prog,
        random_output_prog
#endif
    };

    // Compile and engine creation
    auto executable = poplar::compileGraph(
        graph, vector<poplar::program::Program>{n_random});

    // Create an engine from the executable.
    poplar::Engine engine{move(executable)};

    // Connect I/O to the engine.
    engine.connectStreamToCallback("random_output_stream",
        unique_ptr<random_output_stream_callback>{new random_output_stream_callback()});

    // Attach to first available device.
    poplar::Device device;
    bool attached = false;
    for (auto& device_candidate : devices) {
        if (device_candidate.attach()) {
            device = move(device_candidate);
            attached = true;
            break;
        }
    }

    if (!attached) {
        throw runtime_error("Failed to attach to a device");
    }

    // Load the engine to the device.
    engine.load(device);

    // Run it!
    engine.run(0);
}

int main(int argc, char*argv[]) {
    (void)argc;
    (void)argv;

    if (argc == 2) {
      use_random = true;
      random_seed =  atoi(argv[1]);
      std::cout << "user arg : seed state " << random_seed << std::endl;
    }
    else if (argc == 4) {
      use_random = false;
      base = atoi(argv[1]);
      tile_factor = atoi(argv[2]);
      worker_factor = atoi(argv[3]);
      std::cout << "user arg : base " << base << ", tile factor " << tile_factor << ",worker_factor " << worker_factor << std::endl;
    }

    if (use_random)
      std::cout << "using random seed state " << random_seed << std::endl;
    else
      std::cout << "using non-random seed state  base " << base << " tile factor " << tile_factor << " worker_factor " << worker_factor << std::endl;

    ipu_program(random_output_floats);

    return 0;
}
