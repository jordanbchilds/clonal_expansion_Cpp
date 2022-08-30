#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poplar/Graph.hpp>
#include "codelets/random.hpp"

using namespace std;

const size_t per_worker = 8;
const size_t workers = 6;
const size_t random_output_floats = workers * per_worker;

class random_output_stream_callback : public poplar::StreamCallback {
    public:
        random_output_stream_callback(void* pData, uint32_t bytes) :
            pData(static_cast<uint8_t*>(pData)), bytes(bytes) {}

        void fetch(void *p) override {
            // Just print results.
            float* random_output = static_cast<float*>(p);
            cout << endl << "Results:" << endl;
            for (uint32_t e = 0; e < random_output_floats;) {
              for (uint32_t r = 0; r <  per_worker; ++r) {
                 cout << random_output[e] << " ";
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

        uint8_t* pData;
        uint32_t bytes;
};

void ipu_program(float* random_output_ptr, size_t random_output_floats) {
    const uint32_t num_ipus = 1;
    const uint32_t num_tiles = 1;

    // Enumerate IPUs.
    poplar::DeviceManager manager{};
    auto devices = manager.getDevices(poplar::TargetType::IPU, num_ipus);
    if (devices.empty()) {
        throw poplar::runtime_error("No IPU devices available num_ipus=" + num_ipus);
    }

    auto& target = devices[0].getTarget();

    //auto graph = unique_ptr<poplar::Graph>(new poplar::Graph{target});
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

    // Program part for compute.
    auto random_prog = poplar::program::Sequence{};

    // Add op to program.
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
        random_prog,
        random_output_prog,
        random_prog,
        random_output_prog,
        random_prog,
        random_output_prog
    };

    // Compile and engine creation
    auto executable = poplar::compileGraph(
        graph, vector<poplar::program::Program>{n_random});

    // Create an engine from the executable.
    poplar::Engine engine{move(executable)};

    // Connect I/O to the engine.
    engine.connectStreamToCallback("random_output_stream",
        unique_ptr<random_output_stream_callback>{new random_output_stream_callback(
            random_output_ptr, random_output_floats * sizeof(poplar::FLOAT))});

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
    vector<float> random_output(random_output_floats, 0.0f);
    ipu_program(random_output.data(), random_output.size());
    return 0;
}
