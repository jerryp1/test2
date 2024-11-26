#include "BuildingBlocks/geometric_perspective_protocols.h"
#include "BuildingBlocks/value-extension.h"
#include "BuildingBlocks/truncation.h"
#include <iostream>

using namespace sci;
using namespace std;

int dim, bwA, bwB, shift;
// vars
int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
GeometricPerspectiveProtocols *geom;
Truncation *trunc_oracle;
AuxProtocols *aux;
XTProtocol *xtProtocol;
PRG128 prg;

void sirnn_unsigned_mul() {
    uint64_t *inA = new uint64_t[dim];
    uint64_t *inB = new uint64_t[dim];
    uint64_t *outC = new uint64_t[dim];
    prg.random_data(inA, dim * sizeof(uint64_t));
    prg.random_data(inB, dim * sizeof(uint64_t));
    uint64_t mask_bwA = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    uint64_t mask_bwB = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
    for (int i = 0; i < dim; i++) {
        inA[i] = inA[i] & mask_bwA;
        inB[i] = inB[i] & mask_bwB;
        outC[i] = 0;
    }
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    geom->sirnn_unsigned_mul(dim, inA, inB, outC, bwA, bwB, bwA + bwB);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for SirNN mul = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "unsigned mul " << bwA << " " << bwB << " bits ints. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        uint64_t mask_out = (1ULL << (bwA + bwB)) - 1;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *inB_bob = new uint64_t[dim];
        uint64_t *res = new uint64_t[dim];
        uint64_t *outC_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(inB_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outC_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bwA;
            inB[i] = (inB[i] + inB_bob[i]) & mask_bwB;
            res[i] = (inA[i] * inB[i]) & mask_out;
            outC[i] = (outC[i] + outC_bob[i]) & mask_out;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            assert(res[i] == outC[i]);
            if (res[i] == outC[i]) {
                count++;
            }
        }
        cout << count << " / " << dim << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(inB, sizeof(uint64_t) * dim);
        iopack->io->send_data(outC, sizeof(uint64_t) * dim);
    }
}

void geom_signed_mul() {
    uint64_t *inA = new uint64_t[dim];
    uint64_t *inB = new uint64_t[dim];
    uint64_t mask_bwA = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    uint64_t mask_bwB = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
    if (party == sci::ALICE) {
        uint64_t *inputA = new uint64_t[dim];
        uint64_t *inputB = new uint64_t[dim];
        uint64_t maskA = (1ULL << (bwA - 2)) - 1;
        uint64_t maskB = (1ULL << (bwB - 2)) - 1;
        prg.random_data(inputA, dim * sizeof(uint64_t));
        prg.random_data(inputB, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inputA[i] = inputA[i] & maskA;
            inputB[i] = inputB[i] & maskB;
        }
        int size1 = dim / 4;
        int size2 = dim / 2;
        int size3 = dim * 3 / 4;
        for (int i = 0; i < size1; i++) {
            inputA[i] = inputA[i];
            inputB[i] = inputB[i];
        }
        for (int i = size1; i < size2; i++) {
            inputA[i] = inputA[i];
            inputB[i] = (1ULL << bwB) - inputB[i];
        }
        for (int i = size2; i < size3; i++) {
            inputA[i] = (1ULL << bwA) - inputA[i];
            inputB[i] = inputB[i];
        }
        for (int i = size3; i < dim; i++) {
            inputA[i] = (1ULL << bwA) - inputA[i];
            inputB[i] = (1ULL << bwB) - inputB[i];
        }
        prg.random_data(inA, dim * sizeof(uint64_t));
        prg.random_data(inB, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inA[i] = inA[i] & mask_bwA;
            inB[i] = inB[i] & mask_bwB;
        }
        uint64_t *bob_inputA = new uint64_t[dim];
        uint64_t *bob_inputB = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            bob_inputA[i] = (inputA[i] - inA[i]) & mask_bwA;
            bob_inputB[i] = (inputB[i] - inB[i]) & mask_bwB;
        }
        iopack->io->send_data(bob_inputA, sizeof(uint64_t) * dim);
        iopack->io->send_data(bob_inputB, sizeof(uint64_t) * dim);
    } else {
        iopack->io->recv_data(inA, sizeof(uint64_t) * dim);
        iopack->io->recv_data(inB, sizeof(uint64_t) * dim);
    }
    uint64_t *outC = new uint64_t[dim];
    uint64_t mask_out = (1ULL << (bwA + bwB)) - 1;
    for (int i = 0; i < dim; i++) {
        outC[i] = 0;
    }
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    geom->signed_mul(dim, inA, inB, outC, bwA, bwB, bwA + bwB);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Geometric perspective mul = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "Geometric Perspective signed mul " << bwA << " " << bwB << " bits ints. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *inB_bob = new uint64_t[dim];
        uint64_t *res = new uint64_t[dim];
        uint64_t *outC_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(inB_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outC_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bwA;
            inB[i] = (inB[i] + inB_bob[i]) & mask_bwB;
            res[i] = (signed_val(inA[i], bwA) * signed_val(inB[i], bwB)) & mask_out;
            outC[i] = (outC[i] + outC_bob[i]) & mask_out;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            assert(res[i] == outC[i]);
            if (res[i] == outC[i]) {
                count++;
            }
        }
        cout << count << " / " << dim << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(inB, sizeof(uint64_t) * dim);
        iopack->io->send_data(outC, sizeof(uint64_t) * dim);
    }
}

void geom_sign_extension() {
    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    if (party == sci::ALICE) {
        uint64_t *inputA = new uint64_t[dim];
        uint64_t mask = (1ULL << (bwA - 2)) - 1;
        prg.random_data(inputA, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inputA[i] = inputA[i] & mask;
        }
        int size = dim / 2;
        for (int i = size; i < dim; i++) {
          inputA[i] = (1ULL << bwA) - inputA[i];
        }
        prg.random_data(inA, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inA[i] &= mask_bw;
            outB[i] = 0;
        }
        uint64_t *bob_input = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            bob_input[i] = (inputA[i] - inA[i]) & mask_bw;
        }
        iopack->io->send_data(bob_input, sizeof(uint64_t) * dim);
    } else {
        iopack->io->recv_data(inA, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            outB[i] = 0;
        }
    }
    uint64_t mask_out = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    geom->sign_extension(dim, inA, outB, bwA, bwB);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Geometric Perspective value extension = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "Geometric Perspective sign extension " << bwA << " " << bwB << " bits ints bits. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *outB_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bw;
            outB[i] = (outB[i] + outB_bob[i]) & mask_out;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            if (signed_val(outB[i], bwB) == signed_val(inA[i], bwA)) {
                count++;
            }
        }
        cout << count << " / " << dim << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(outB, sizeof(uint64_t) * dim);
    }
}

void sirnn_sign_extension() {
    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    prg.random_data(inA, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
        inA[i] &= mask_bw;
        outB[i] = 0;
    }
    uint64_t mask_out = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
    int64_t c0 = xtProtocol->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    xtProtocol->s_extend(dim, inA, outB, bwA, bwB);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for SirNN value extension = " << (temp / 1000.0) << std::endl;
    int64_t c1 = xtProtocol->iopack->io->counter;
    std::cout << "SirNN sign extension " << bwA << " " << bwB << " bits ints bits. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *outB_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bw;
            outB[i] = (outB[i] + outB_bob[i]) & mask_out;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            if (signed_val(outB[i], bwB) == signed_val(inA[i], bwA)) {
                count++;
            }
        }
        cout << count << " / " << dim << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(outB, sizeof(uint64_t) * dim);
    }
}

void geom_trunc() {
    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    if (party == sci::ALICE) {
        uint64_t *inputA = new uint64_t[dim];
        uint64_t mask = (1ULL << (bwA - 2)) - 1;
        prg.random_data(inputA, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inputA[i] = inputA[i] & mask;
        }
        int size = dim / 2;
        for (int i = size; i < dim; i++) {
            inputA[i] = (1ULL << bwA) - inputA[i];
        }
        prg.random_data(inA, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inA[i] &= mask_bw;
            outB[i] = 0;
        }
        uint64_t *bob_input = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            bob_input[i] = (inputA[i] - inA[i]) & mask_bw;
        }
        iopack->io->send_data(bob_input, sizeof(uint64_t) * dim);
    } else {
        iopack->io->recv_data(inA, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            outB[i] = 0;
        }
    }
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    geom->new_truncate(dim, inA, outB, shift, bwA);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Geometric Perspective Protocol = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "truncate " << "signed" << " " << bwA << " bits ints ";
    std::cout <<  "shift by " << shift << " bits. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *outB_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bw;
            outB[i] = (outB[i] + outB_bob[i]) & mask_bw;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            assert((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA));
            if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA)) {
                count++;
            }
        }
        cout << count << " / " << dim << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(outB, sizeof(uint64_t) * dim);
    }
}

void trunc(bool signed_arithmetic = true) {
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    prg.random_data(inA, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
        inA[i] &= mask_bw;
        outB[i] = 0;
    }
    int64_t c0 = trunc_oracle->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    trunc_oracle->truncate(dim, inA, outB, shift, bwA, signed_arithmetic);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for SirNN = " << (temp / 1000.0) << std::endl;
    int64_t c1 = trunc_oracle->iopack->io->counter;
    std::cout << "truncate " << (signed_arithmetic ? "signed" : "nonsigned") << " " << bwA << " bits ints ";
    std::cout <<  "shift by " << shift << " bits. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *outB_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bw;
            outB[i] = (outB[i] + outB_bob[i]) & mask_bw;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            if (signed_arithmetic) {
                if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA)) {
                    count++;
                }
            } else {
                if ((inA[i] >> shift) == outB[i]) {
                    count++;
                }
            }
        }
        cout << count << " / " << dim << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(outB, sizeof(uint64_t) * dim);
    }
}

void msb0_trunc() {
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    prg.random_data(inA, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
        inA[i] &= mask_bw;
        outB[i] = 0;
    }
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    geom->msb0_truncation(dim, inA, outB, shift, bwA);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Cheetah = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "truncate " << " " << bwA << " bits ints ";
    std::cout <<  "shift by " << shift << " bits. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        int total = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *outB_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bw;
            outB[i] = (outB[i] + outB_bob[i]) & mask_bw;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            if (inA[i] < (1ULL << (bwA - 2))) {
                total++;
                if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA)) {
                    count++;
                }
                if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA) + 1) {
                    count++;
                }
                if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA) - 1) {
                    count++;
                }
            }
        }
        cout << count << " / " << total << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(outB, sizeof(uint64_t) * dim);
    }
}

void geom_trunc_with_one_bit_error() {
    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    if (party == sci::ALICE) {
        uint64_t *inputA = new uint64_t[dim];
        uint64_t mask = (1ULL << (bwA - 2)) - 1;
        prg.random_data(inputA, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inputA[i] = inputA[i] & mask;
        }
        int size = dim / 2;
        for (int i = size; i < dim; i++) {
            inputA[i] = (1ULL << bwA) - inputA[i];
        }
        prg.random_data(inA, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            inA[i] &= mask_bw;
            outB[i] = 0;
        }
        uint64_t *bob_input = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            bob_input[i] = (inputA[i] - inA[i]) & mask_bw;
        }
        iopack->io->send_data(bob_input, sizeof(uint64_t) * dim);
    } else {
        iopack->io->recv_data(inA, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            outB[i] = 0;
        }
    }
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    geom->truncate_with_one_bit_error(dim, inA, outB, shift, bwA);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Geometric Perspective truncation with one bit error = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "truncate " << " " << bwA << " bits ints ";
    std::cout <<  "shift by " << shift << " bits. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        int total = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *outB_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bw;
            outB[i] = (outB[i] + outB_bob[i]) & mask_bw;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA)) {
                count++;
            }
            if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA) + 1) {
                count++;
            }
            if ((signed_val(inA[i], bwA) >> shift) == signed_val(outB[i], bwA) - 1) {
                count++;
            }
        }
        cout << count << " / " << dim << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(outB, sizeof(uint64_t) * dim);
    }
}

void geom_msb0_trunc() {
    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    prg.random_data(inA, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
        inA[i] &= mask_bw;
        outB[i] = 0;
    }
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    geom->msb0_truncation(dim, inA, outB, shift, bwA);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Geometric Perspective truncation with msb0 = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "truncate " << " " << bwA << " bits ints ";
    std::cout <<  "shift by " << shift << " bits. Sent " << ((c1 - c0) * 8 / dim) << " bits\n";
    if (party == sci::ALICE) {
        int count = 0;
        int total = 0;
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t *outB_bob = new uint64_t[dim];
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < dim; i++) {
            inA[i] = (inA[i] + inA_bob[i]) & mask_bw;
            outB[i] = (outB[i] + outB_bob[i]) & mask_bw;
            cout << inA[i] << " "  << (inA[i] >> shift) << " " << outB[i] << endl;
        }
        cout << "Testing for correctness..." << endl;
        for (int i = 0; i < dim; i++) {
            if (inA[i] < (1ULL << (bwA -1))) {
                total++;
                if (outB[i] == (inA[i] >> shift)) {
                    count++;
                }
            }
        }
        cout << count << " / " << total << " Correct!" << endl;
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        iopack->io->send_data(outB, sizeof(uint64_t) * dim);
    }
}

void test_euclidean_distance(ReLUProtocol<uint64_t> *maxpool) {
    int size = 10;
    uint64_t *inA = new uint64_t[dim];
    uint64_t **inB = new uint64_t*[size];
    for (int i = 0; i < size; i++) {
      inB[i] = new uint64_t[dim];
    }
    uint64_t mask_bw = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    uint64_t maskA = (1ULL << (bwA - 2)) - 1;

    if (party == sci::ALICE) {
        uint64_t *inputA = new uint64_t[dim];
        uint64_t **inputB = new uint64_t*[size];
        for (int i = 0; i < size; i++) {
          inputB[i] = new uint64_t[dim];
        }
        prg.random_data(inputA, dim * sizeof(uint64_t));
        for (int i = 0; i < size; i++) {
          prg.random_data(inputB[i], dim * sizeof(uint64_t));
        }
        for (int i = 0; i < dim; i++) {
            inputA[i] = inputA[i] & maskA;
        }
        for (int i = 0; i < size; i++) {
          for (int j = 0; j < dim; j++) {
            inputB[i][j] = inputB[i][j] & maskA;
          }
        }
        int size1 = dim / 4;
        int size2 = dim / 2;
        int size3 = dim * 3 / 4;
        for (int i = 0; i < size1; i++) {
            inputA[i] = inputA[i];
        }
        for (int i = 0; i < size; i++) {
          for (int j = 0; j < size1; j++) {
            inputB[i][j] = inputB[i][j];
          }
        }
        for (int i = size1; i < size2; i++) {
            inputA[i] = inputA[i];
        }
        for (int i = 0; i < size; i++) {
            for (int j = size1; j < size2; j++) {
                inputB[i][j] = (1ULL << bwB) - inputB[i][j];
            }
        }
        for (int i = size2; i < size3; i++) {
            inputA[i] = (1ULL << bwA) - inputA[i];
        }
        for (int i = 0; i < size; i++) {
            for (int j = size2; j < size3; j++) {
                inputB[i][j] = inputB[i][j];
            }
        }
        for (int i = size3; i < dim; i++) {
            inputA[i] = (1ULL << bwA) - inputA[i];
        }
        for (int i = 0; i < size; i++) {
            for (int j = size3; j < dim; j++) {
                inputB[i][j] = (1ULL << bwB) - inputB[i][j];
            }
        }
        prg.random_data(inA, dim * sizeof(uint64_t));
        for (int i = 0; i < size; i++) {
          prg.random_data(inB[i], dim * sizeof(uint64_t));
        }
        for (int i = 0; i < dim; i++) {
            inA[i] = inA[i] & maskA;
        }
        for (int i = 0; i < size; i++) {
          for (int j = 0; j < dim; j++) {
            inB[i][j] = inB[i][j] & maskA;
          }
        }
        uint64_t *bob_inputA = new uint64_t[dim];
        uint64_t **bob_inputB = new uint64_t*[size];
        for (int i = 0; i < size; i++) {
          bob_inputB[i] = new uint64_t[dim];
        }
        for (int i = 0; i < size; i++) {
          for (int j = 0; j < dim; j++) {
              bob_inputB[i][j] = (inputB[i][j] - inB[i][j]) & maskA;
          }
        }
        for (int i = 0; i < dim; i++) {
            bob_inputA[i] = (inputA[i] - inA[i]) & maskA;
        }
        iopack->io->send_data(bob_inputA, sizeof(uint64_t) * dim);
        for (int i = 0; i < size; i++) {
            iopack->io->send_data(bob_inputB[i], sizeof(uint64_t) * dim);
        }
    } else {
        iopack->io->recv_data(inA, sizeof(uint64_t) * dim);
        for (int i = 0; i < size; i++) {
            iopack->io->recv_data(inB[i], sizeof(uint64_t) * dim);
        }
    }
    uint64_t mask_out = (1ULL << (bwA + bwA)) - 1;
    int64_t c0 = geom->iopack->io->counter;
    std::cout << c0 << std::endl;
    auto start_timer = std::chrono::high_resolution_clock::now();
    uint64_t outC = geom->euclidean_distance(size, dim, inA, inB, bwA, maxpool);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Euclidean distance = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << c1 << std::endl;
    std::cout << "Euclidean distance test " << bwA << " bits ints. Sent " << ((c1 - c0) * 8) << " bits\n";
    if (party == sci::ALICE) {
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t **inB_bob = new uint64_t*[size];
        for (int i = 0; i < size; i++) {
          inB_bob[i] = new uint64_t[dim];
        }
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < size; i++) {
            iopack->io->recv_data(inB_bob[i], sizeof(uint64_t) * dim);
        }
        uint64_t *outC_bob = new uint64_t[1];
        iopack->io->recv_data(outC_bob, sizeof(uint64_t) * 1);
        outC = (outC_bob[0] + outC) & mask_out;
        uint64_t res = -1;
        for (int i = 0; i < size; i++) {
            uint64_t temp = 0;
            for (int j = 0; j < dim; j++) {
                uint64_t a = (inA[j] + inA_bob[j]) & mask_bw;
                uint64_t b = (inB[i][j] + inB_bob[i][j]) & mask_bw;
                temp = (temp + (a * b)) & mask_out;
            }
            if (temp < res) {
              res = temp;
            }
        }
        cout << "Testing for correctness..." << endl;
        //assert(res == outC);
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        for (int i = 0; i < size; i++) {
            iopack->io->send_data(inB[i], sizeof(uint64_t) * dim);
        }
        uint64_t *inA_bob = new uint64_t[1];
        inA_bob[0] = outC;
        iopack->io->send_data(inA_bob, sizeof(uint64_t) * 1);
    }
}

void test_sirnn_euclidean_distance(ReLUProtocol<uint64_t> *maxpool) {
    int size = 10;
    uint64_t *inA = new uint64_t[dim];
    uint64_t **inB = new uint64_t*[size];
    for (int i = 0; i < size; i++) {
        inB[i] = new uint64_t[dim];
    }
    prg.random_data(inA, dim * sizeof(uint64_t));
    for (int i = 0; i < size; i++) {
        prg.random_data(inB[i], dim * sizeof(uint64_t));
    }
    uint64_t mask_bwA = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    for (int i = 0; i < dim; i++) {
        inA[i] = inA[i] & mask_bwA;
    }
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < dim; j++) {
        inB[i][j] = inB[i][j] & mask_bwA;
      }
    }
    uint64_t mask_out = (1ULL << (bwA + bwA)) - 1;
    int64_t c0 = geom->iopack->io->counter;
    auto start_timer = std::chrono::high_resolution_clock::now();
    uint64_t outC = geom->sirnn_euclidean_distance(size, dim, inA, inB, bwA, maxpool);
    auto temp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for Euclidean distance = " << (temp / 1000.0) << std::endl;
    int64_t c1 = geom->iopack->io->counter;
    std::cout << "Euclidean distance test " << bwA << " bits ints. Sent " << ((c1 - c0) * 8) << " bits\n";
    if (party == sci::ALICE) {
        uint64_t *inA_bob = new uint64_t[dim];
        uint64_t **inB_bob = new uint64_t*[size];
        for (int i = 0; i < size; i++) {
          inB_bob[i] = new uint64_t[dim];
        }
        iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
        for (int i = 0; i < size; i++) {
            iopack->io->recv_data(inB_bob[i], sizeof(uint64_t) * dim);
        }
        uint64_t *outC_bob = new uint64_t[1];
        iopack->io->recv_data(outC_bob, sizeof(uint64_t) * 1);
        outC = (outC_bob[0] + outC) & mask_out;
        uint64_t res = -1;
        for (int i = 0; i < size; i++) {
            uint64_t temp = 0;
            for (int j = 0; j < dim; j++) {
                uint64_t a = (inA[j] + inA_bob[j]) & mask_bwA;
                uint64_t b = (inB[i][j] + inB_bob[i][j]) & mask_bwA;
                temp = (temp + (a * b)) & mask_out;
            }
            if (temp < res) {
              res = temp;
            }
        }
        cout << "Testing for correctness..." << endl;
        //assert(res == outC);
    } else { // BOB
        iopack->io->send_data(inA, sizeof(uint64_t) * dim);
        for (int i = 0; i < size; i++) {
            iopack->io->send_data(inB[i], sizeof(uint64_t) * dim);
        }
        uint64_t *inA_bob = new uint64_t[1];
        inA_bob[0] = outC;
        iopack->io->send_data(inA_bob, sizeof(uint64_t) * 1);
    }
}

int main(int argc, char **argv) {
    // default
    dim = 1 << 16;
    bwA = 25;
    bwB = 30;
    shift = 12;

    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of operations");
    amap.arg("bwA", bwA, "Bitlength of inputA");
    amap.arg("bwB", bwB, "Bitlength of inputB");
    amap.arg("s", shift, "Bitlength of shift");
    amap.arg("ip", address, "IP Address of server (ALICE)");

    amap.parse(argc, argv);

    auto start_timer = std::chrono::high_resolution_clock::now();

    iopack = new IOPack(party, port, "127.0.0.1");
    otpack = new OTPack(iopack, party);
    geom = new GeometricPerspectiveProtocols(party, iopack, otpack);
    trunc_oracle = new Truncation(party, iopack, otpack);
    aux = new AuxProtocols(party, iopack, otpack);
    xtProtocol = new XTProtocol(party, iopack, otpack);
    auto temp = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start_timer).count();
    std::cout << "Time in sec for preprocessing = " << (temp / 1000.0) << std::endl;

//    cout << "<><><><> test geometric perspective faithful truncation <><><><>" << endl;
//    geom_trunc();
//
//    cout << "<><><><> test geometric perspective signed multiplication <><><><>" << endl;
//    geom_signed_mul();
//
//    cout << "<><><><> test geometric perspective value extension <><><><>" << endl;
//    geom_sign_extension();
//
//    cout << "<><><><> test geometric perspective truncation with one bit error <><><><>" << endl;
//    geom_trunc_with_one_bit_error();


      geom_msb0_trunc();
//    ReLURingProtocol<uint64_t> *relu = new ReLURingProtocol<uint64_t>(party, RING, iopack, bwA + bwA, MILL_PARAM, otpack);
//    test_euclidean_distance(relu);
//
//    test_sirnn_euclidean_distance(relu);
}
