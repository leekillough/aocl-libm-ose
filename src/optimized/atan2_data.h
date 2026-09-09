/*
 * Copyright (C) 2026 Advanced Micro Devices, Inc. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef ATAN2_DATA_H
#define ATAN2_DATA_H

#ifdef ATAN2_VRD
/* Shared constants for vrd2/vrd4/vrd8 atan2 kernels. */
/*
 * Scaling to handle overflow / underflow :
 * The reduced-region terms are pre-scaled by 1/4. rr is a ratio, so an exact
 * power-of-two scale leaves it bit-identical, but it keeps den*sqrt3 + num
 * from overflowing once max(|x|,|y|) passes DBL_MAX/(1+sqrt3) ~ 6.58e307.
 * SMALL_LIM is the matching lower guard: below it the operands are scaled up
 * first, otherwise the 1/4 would push the reduced-region terms subnormal.
 */
#define ATAN2_VRD_SMALL_LIM   0x0050000000000000UL  /* 2^-1018     */
#define ATAN2_VRD_PRESCALE    0x1p-2                /* 1/4         */
#define ATAN2_VRD_SQRT3_PS    0x1.bb67ae8584caap-2  /* sqrt(3)/4   */

#define ATAN2_VRD_RANGE       0x1.126145e9ecd56p-2  /* 2 - sqrt(3) */
#define ATAN2_VRD_PI6         0x1.0c152382d7366p-1  /* pi/6        */
#define ATAN2_VRD_PI2         0x1.921fb54442d18p0   /* pi/2        */
#define ATAN2_VRD_PI          0x1.921fb54442d18p1   /* pi          */
#define ATAN2_VRD_SCALE_UP    0x1p54
/* atan(rr) = rr + rr^3*(C0 + C1*rr^2 + ...) : coeffs of rr^3..rr^17 */
#define ATAN2_VRD_POLY_C0    (-0x1.55555555553a5p-2)
#define ATAN2_VRD_POLY_C1     0x1.999999993b043p-3
#define ATAN2_VRD_POLY_C2    (-0x1.249248eb994b7p-3)
#define ATAN2_VRD_POLY_C3     0x1.c71c506f7a87ep-4
#define ATAN2_VRD_POLY_C4    (-0x1.7457ab2dfd3d2p-4)
#define ATAN2_VRD_POLY_C5     0x1.3a8f4b097ef8dp-4
#define ATAN2_VRD_POLY_C6    (-0x1.09a744e6dd6ap-4)
#define ATAN2_VRD_POLY_C7     0x1.6c4e437d4b353p-5
#endif /* ATAN2_VRD */

#ifdef ATAN2_JBY256
  /* Arrays atan_jby256_lead and atan_jby256_tail contain
     leading and trailing parts respectively of precomputed
     values of atan(j/256), for j = 16, 17, ..., 256.
     atan_jby256_lead contains the first 21 bits of precision,
     and atan_jby256_tail contains a further 53 bits precision. */
static struct {
    double head;
    double tail;
} atan_jby256[241] = {

{0x1.ff55bp-5, 0x1.6e59fbd38db2dp-26},
{0x1.0f99ep-4, 0x1.4e3aa54dedf96p-25},
{0x1.1f86dp-4, 0x1.7e105ab1bda89p-25},
{0x1.2f719p-4, 0x1.8c5254d013fd1p-27},
{0x1.3f59fp-4, 0x1.cf8ab3ad62671p-29},
{0x1.4f3fdp-4, 0x1.9dca4bec8046ap-26},
{0x1.5f232p-4, 0x1.3f4b5ec98a8dap-26},
{0x1.6f03bp-4, 0x1.b9d49619d81fep-25},
{0x1.7ee18p-4, 0x1.3017887460935p-27},
{0x1.8ebc5p-4, 0x1.11e3eca0b9944p-26},
{0x1.9e941p-4, 0x1.4f3f73c5a332ep-26},
{0x1.ae68ap-4, 0x1.c71c8ae0e00a7p-26},
{0x1.be39ep-4, 0x1.7cde0f86fbdc8p-25},
{0x1.ce07cp-4, 0x1.70f328c889c73p-26},
{0x1.ddd21p-4, 0x1.c07ae9b994fp-26},
{0x1.ed98cp-4, 0x1.0c8021d7b169ap-27},
{0x1.fd5bap-4, 0x1.35585edb8cb22p-25},
{0x1.068d5p-3, 0x1.0842567b30e97p-24},
{0x1.0e6adp-3, 0x1.99e811031472ep-24},
{0x1.16465p-3, 0x1.041821416bcefp-25},
{0x1.1e1fap-3, 0x1.f6086e4dc96f5p-24},
{0x1.25f6ep-3, 0x1.71a535c5f1b5fp-27},
{0x1.2dcbdp-3, 0x1.65f743fe63ca2p-24},
{0x1.359e8p-3, 0x1.dbd733472d014p-24},
{0x1.3d6eep-3, 0x1.d18cc4d8b0d1ep-24},
{0x1.453cep-3, 0x1.8c12553c8fb2ap-24},
{0x1.4d087p-3, 0x1.53b49e2e8f992p-24},
{0x1.54d18p-3, 0x1.7422ae148c141p-24},
{0x1.5c981p-3, 0x1.e3ec269df56aap-27},
{0x1.645bfp-3, 0x1.ff6754e7e0ac9p-24},
{0x1.6c1d4p-3, 0x1.131267b1b5aadp-24},
{0x1.73dbdp-3, 0x1.d14fa403a94bcp-24},
{0x1.7b97bp-3, 0x1.2f396c089a3d8p-25},
{0x1.8350bp-3, 0x1.c731d78fa95bbp-24},
{0x1.8b06ep-3, 0x1.c50f385177399p-24},
{0x1.92ba3p-3, 0x1.f41409c6f2c21p-25},
{0x1.9a6a8p-3, 0x1.d2d90c4c39ecp-24},
{0x1.a217ep-3, 0x1.80420696f2107p-25},
{0x1.a9c23p-3, 0x1.b40327943a2fp-27},
{0x1.b1696p-3, 0x1.5d35e02f3d2a4p-25},
{0x1.b90d7p-3, 0x1.4a498288117b1p-25},
{0x1.c0ae5p-3, 0x1.35da119afb324p-25},
{0x1.c84bfp-3, 0x1.14e85cdb9a909p-24},
{0x1.cfe65p-3, 0x1.38754e5547b9ap-25},
{0x1.d77d5p-3, 0x1.be40ae6ce3247p-24},
{0x1.df11p-3, 0x1.0c993b3bea7e7p-24},
{0x1.e6a14p-3, 0x1.1d2dd89ac3359p-24},
{0x1.ee2e1p-3, 0x1.1476603332c48p-25},
{0x1.f5b75p-3, 0x1.f25901bac55b7p-24},
{0x1.fd3d1p-3, 0x1.f881b7c826e29p-24},
{0x1.025fap-2, 0x1.441996d698d21p-24},
{0x1.061eep-2, 0x1.407ac521ea08ap-23},
{0x1.09dc5p-2, 0x1.2fb0c6c4b1724p-23},
{0x1.0d97ep-2, 0x1.ca135966a3e19p-23},
{0x1.1151ap-2, 0x1.b1218e4d646e4p-25},
{0x1.15097p-2, 0x1.d4e72a350d28bp-25},
{0x1.18bf5p-2, 0x1.4617e2f04c329p-23},
{0x1.1c735p-2, 0x1.096ec41e82654p-25},
{0x1.20255p-2, 0x1.9f91f25773e6ep-24},
{0x1.23d56p-2, 0x1.59c0820f1d674p-25},
{0x1.27837p-2, 0x1.02bf7a2df1064p-25},
{0x1.2b2f7p-2, 0x1.fb36bfc40508fp-23},
{0x1.2ed98p-2, 0x1.ea08f3f8dc893p-24},
{0x1.32818p-2, 0x1.3ed6254656a0fp-24},
{0x1.36277p-2, 0x1.b83f5e5e69c5bp-25},
{0x1.39cb4p-2, 0x1.d6ec2af768593p-23},
{0x1.3d6d1p-2, 0x1.493889a226f95p-25},
{0x1.410cbp-2, 0x1.5ad8fa65279bap-23},
{0x1.44aa4p-2, 0x1.b615784d45438p-25},
{0x1.4845ap-2, 0x1.09a184368f145p-23},
{0x1.4bdeep-2, 0x1.61a2439b0d91cp-24},
{0x1.4f75fp-2, 0x1.ce1a65e39a979p-24},
{0x1.530adp-2, 0x1.32a39a93b6a67p-23},
{0x1.569d8p-2, 0x1.1c3699af804e8p-23},
{0x1.5a2ep-2, 0x1.75e0f4e44edeap-26},
{0x1.5dbc3p-2, 0x1.f77ced1a7a83bp-23},
{0x1.61484p-2, 0x1.84e7f0cb1b51dp-29},
{0x1.64d1fp-2, 0x1.ec6b838b02dfep-23},
{0x1.68597p-2, 0x1.3ebf4dfbeda88p-23},
{0x1.6bdeap-2, 0x1.9397aed9cb475p-23},
{0x1.6f619p-2, 0x1.07937bc239c55p-24},
{0x1.72e22p-2, 0x1.aa754553131b6p-23},
{0x1.76607p-2, 0x1.4a05d407c45ddp-24},
{0x1.79dc6p-2, 0x1.132231a206ddp-23},
{0x1.7d56p-2, 0x1.2d8ecfdd69c88p-24},
{0x1.80cd4p-2, 0x1.a852c74218606p-24},
{0x1.84422p-2, 0x1.71bf2baeebb5p-23},
{0x1.87b4bp-2, 0x1.83d7db749182bp-27},
{0x1.8b24dp-2, 0x1.ca50d92b6da17p-25},
{0x1.8e929p-2, 0x1.6f5cde85302ap-26},
{0x1.91fdep-2, 0x1.f343198910742p-24},
{0x1.9566dp-2, 0x1.0e8d241ccd80bp-24},
{0x1.98cd5p-2, 0x1.1535ac619e6cap-24},
{0x1.9c316p-2, 0x1.7316041c36cd3p-24},
{0x1.9f93p-2, 0x1.985a000637d9p-24},
{0x1.a2f23p-2, 0x1.f2f29858c0a6bp-25},
{0x1.a64eep-2, 0x1.879847f96d90ap-23},
{0x1.a9a92p-2, 0x1.ab3d319e12e43p-23},
{0x1.ad00fp-2, 0x1.5088162dfc4c3p-24},
{0x1.b0564p-2, 0x1.05749a1cd9d8dp-25},
{0x1.b3a91p-2, 0x1.da65c6c6b861fp-26},
{0x1.b6f96p-2, 0x1.739bf7df1ad66p-25},
{0x1.ba473p-2, 0x1.bc31252aa3343p-25},
{0x1.bd928p-2, 0x1.e528191ad3aafp-26},
{0x1.c0db4p-2, 0x1.929d93df19f19p-23},
{0x1.c4219p-2, 0x1.ff11eb693a083p-26},
{0x1.c7655p-2, 0x1.55ae3f145a3ap-27},
{0x1.caa68p-2, 0x1.cbcd8c6c0ca83p-24},
{0x1.cde53p-2, 0x1.0cb04d425d306p-24},
{0x1.d1215p-2, 0x1.9adfcab5be678p-24},
{0x1.d45aep-2, 0x1.93d90c5662508p-23},
{0x1.d791fp-2, 0x1.68489bd35ff41p-24},
{0x1.dac67p-2, 0x1.586ed3da2b7f2p-28},
{0x1.ddf85p-2, 0x1.7604d2e850eefp-23},
{0x1.e127bp-2, 0x1.ac1d12bfb53dap-24},
{0x1.e4548p-2, 0x1.9b3d468274754p-28},
{0x1.e77ebp-2, 0x1.fc5d68d10e53ep-24},
{0x1.eaa65p-2, 0x1.8f9e51884becbp-23},
{0x1.edcb6p-2, 0x1.a87f0869c06d1p-23},
{0x1.f0edep-2, 0x1.31e7279f685fbp-23},
{0x1.f40ddp-2, 0x1.6a8282f9719b5p-27},
{0x1.f72b2p-2, 0x1.0d2724a8a44e1p-25},
{0x1.fa45dp-2, 0x1.a60524b11ad4ep-23},
{0x1.fd5ep-2, 0x1.75fdf832750f6p-26},
{0x1.0039cp-1, 0x1.cf06902e4cd37p-23},
{0x1.01c34p-1, 0x1.e82422d4f6d14p-25},
{0x1.034b7p-1, 0x1.24a091063e6cep-26},
{0x1.04d25p-1, 0x1.8a1a172dc6f39p-24},
{0x1.0657ep-1, 0x1.29b6619f8a92ep-22},
{0x1.07dc3p-1, 0x1.9274d9c1b70cap-24},
{0x1.095f3p-1, 0x1.0c34b1fbb7938p-26},
{0x1.0ae0ep-1, 0x1.639866c20eb58p-25},
{0x1.0c614p-1, 0x1.6d6d0f6832e9fp-23},
{0x1.0de05p-1, 0x1.af54def99f25fp-22},
{0x1.0f5e2p-1, 0x1.16cfc52a00262p-22},
{0x1.10daap-1, 0x1.dcc1e83569c34p-23},
{0x1.1255dp-1, 0x1.37f7a551ed425p-22},
{0x1.13cfbp-1, 0x1.f6360adc98887p-22},
{0x1.15485p-1, 0x1.2c6ec8d35a2c2p-22},
{0x1.16bfap-1, 0x1.bd44df84cb038p-23},
{0x1.1835ap-1, 0x1.117cf826e310ep-22},
{0x1.19aa5p-1, 0x1.ca533f332cfcap-22},
{0x1.1b1dcp-1, 0x1.0f208509dbc2ep-22},
{0x1.1c8fep-1, 0x1.cd07d93c945dep-23},
{0x1.1e00bp-1, 0x1.57bdfd67e6d72p-22},
{0x1.1f704p-1, 0x1.aab89c516c65ap-24},
{0x1.20de8p-1, 0x1.3e823b1a1b8a3p-25},
{0x1.224b7p-1, 0x1.307464a9d6d3ep-23},
{0x1.23b71p-1, 0x1.c5993cd438844p-22},
{0x1.25217p-1, 0x1.ba2fca02ab554p-22},
{0x1.268a9p-1, 0x1.01a5b6983a269p-23},
{0x1.27f26p-1, 0x1.273d1b350efcep-25},
{0x1.2958ep-1, 0x1.64c238c37b0c6p-23},
{0x1.2abe2p-1, 0x1.aded07370a3p-25},
{0x1.2c221p-1, 0x1.78091197eb47fp-23},
{0x1.2d84cp-1, 0x1.4b0f245e0dabep-24},
{0x1.2ee62p-1, 0x1.080d9794e2eafp-22},
{0x1.30464p-1, 0x1.d4ec242b60c77p-23},
{0x1.31a52p-1, 0x1.221d2f940cab8p-27},
{0x1.3302bp-1, 0x1.cdbc42b2bba5ep-24},
{0x1.345fp-1, 0x1.cce37bb440845p-25},
{0x1.35bap-1, 0x1.6c1d999cf1dd1p-22},
{0x1.3713dp-1, 0x1.bed8a07eb0876p-26},
{0x1.386c5p-1, 0x1.69ed88f490e3cp-24},
{0x1.39c39p-1, 0x1.cd41719b73ef3p-25},
{0x1.3b198p-1, 0x1.cbc4ac95b41b7p-22},
{0x1.3c6e4p-1, 0x1.238f1b890f5d8p-22},
{0x1.3dc1cp-1, 0x1.50c4282259cc7p-24},
{0x1.3f13fp-1, 0x1.713d2de87b3e3p-22},
{0x1.4064fp-1, 0x1.1d5a7d2255277p-23},
{0x1.41b4ap-1, 0x1.c0dfd48227ac1p-22},
{0x1.43032p-1, 0x1.1c964dab76753p-22},
{0x1.44506p-1, 0x1.6de56d5704498p-23},
{0x1.459c6p-1, 0x1.4aeb71fd1996ap-23},
{0x1.46e72p-1, 0x1.fbf91c57b1919p-23},
{0x1.4830ap-1, 0x1.d6bef7fbe5d9bp-22},
{0x1.4978fp-1, 0x1.464d3dc249067p-22},
{0x1.4acp-1, 0x1.638e2ec4d9074p-22},
{0x1.4c05ep-1, 0x1.16f4a7247ea7cp-24},
{0x1.4d4a8p-1, 0x1.1a0a740f1d449p-28},
{0x1.4e8dep-1, 0x1.6edbb0114a33ep-23},
{0x1.4fd01p-1, 0x1.dbee8bf1d513fp-24},
{0x1.5111p-1, 0x1.5b8bdb0248f73p-22},
{0x1.5250cp-1, 0x1.7de3d3f5eac65p-22},
{0x1.538f5p-1, 0x1.ee24187ae448bp-23},
{0x1.54ccap-1, 0x1.e06c591ec5192p-22},
{0x1.5608dp-1, 0x1.4e3861a33273ap-24},
{0x1.5743cp-1, 0x1.a9599dcc2bfe6p-24},
{0x1.587d8p-1, 0x1.f732fbad4346cp-25},
{0x1.59b6p-1, 0x1.eb9f573b727dap-22},
{0x1.5aed6p-1, 0x1.8b212a2eb9898p-22},
{0x1.5c239p-1, 0x1.384884c167216p-22},
{0x1.5d589p-1, 0x1.0e2d363020052p-22},
{0x1.5e8c6p-1, 0x1.2820879fbd023p-22},
{0x1.5fbfp-1, 0x1.a1ab9893e4b3p-22},
{0x1.60f08p-1, 0x1.2d1b817a24479p-23},
{0x1.6220dp-1, 0x1.15d7b8ded487bp-25},
{0x1.634ffp-1, 0x1.8968f9db3a5e6p-24},
{0x1.647dep-1, 0x1.71c4171fe136p-22},
{0x1.65aabp-1, 0x1.6d80f605d0d8dp-22},
{0x1.66d66p-1, 0x1.c91f04369159p-24},
{0x1.6800ep-1, 0x1.39f8a15fce2b3p-23},
{0x1.692a4p-1, 0x1.55beda9d94b96p-27},
{0x1.6a527p-1, 0x1.b12c15d60949ap-23},
{0x1.6b798p-1, 0x1.24167b312bfe4p-22},
{0x1.6c9f7p-1, 0x1.0ab8633070278p-22},
{0x1.6dc44p-1, 0x1.54554ebbc80efp-23},
{0x1.6ee7fp-1, 0x1.0204aef5a4bbbp-25},
{0x1.700a7p-1, 0x1.8af08c679cf2dp-22},
{0x1.712bep-1, 0x1.0852a330ae6c8p-22},
{0x1.724c3p-1, 0x1.6d3eb9ec32916p-23},
{0x1.736b6p-1, 0x1.685cb7fcbbaffp-23},
{0x1.74897p-1, 0x1.1f751c1e0bd95p-22},
{0x1.75a67p-1, 0x1.705b1b0f72562p-26},
{0x1.76c24p-1, 0x1.b98d8d808ca93p-22},
{0x1.77dd1p-1, 0x1.2ea22c75cc982p-25},
{0x1.78f6bp-1, 0x1.7aba62bca035p-22},
{0x1.7a0f4p-1, 0x1.d73833442278dp-22},
{0x1.7b26cp-1, 0x1.5a5ca1fb18bf9p-22},
{0x1.7c3d3p-1, 0x1.1a6092b6ecf2dp-25},
{0x1.7d528p-1, 0x1.44fd049aac104p-24},
{0x1.7e66cp-1, 0x1.c114fd8df51ddp-29},
{0x1.7f79ep-1, 0x1.5972f130feae5p-22},
{0x1.808cp-1, 0x1.ca034a55fe19bp-24},
{0x1.819dp-1, 0x1.6e2b149990227p-22},
{0x1.82adp-1, 0x1.b00000294593p-24},
{0x1.83bbep-1, 0x1.8b9bdc442620ep-22},
{0x1.84c9cp-1, 0x1.d94fdfabf3e4fp-23},
{0x1.85d69p-1, 0x1.5db30b145ad9cp-23},
{0x1.86e25p-1, 0x1.e3e1eb95022b1p-23},
{0x1.87edp-1, 0x1.d5b8b45442bd7p-22},
{0x1.88f6bp-1, 0x1.7a046231ecd2ep-22},
{0x1.89ff5p-1, 0x1.feafe3ef55233p-22},
{0x1.8b06fp-1, 0x1.839e7bfd78267p-22},
{0x1.8c0d9p-1, 0x1.45cf49d6fa902p-25},
{0x1.8d132p-1, 0x1.be3132b27f39ap-27},
{0x1.8e17ap-1, 0x1.533980bb84f9fp-22},
{0x1.8f1b3p-1, 0x1.889e2ce3ba395p-26},
{0x1.901dbp-1, 0x1.f7778c3ad0cccp-24},
{0x1.911f3p-1, 0x1.46660cec4eba3p-23},
{0x1.921fbp-1, 0x1.5110b4611a626p-23},
};

#endif /* ATAN2_JBY256 */

typedef struct { double head; double tail; } Atan2LogData;

// 256-entry table: entry (b,k) stores atan(c) as head+tail.
// c(b,k) = 2^(b-4) * (2k+65)/64,  b=0..7, k=0..31
// Reconstructible: asdouble(((uint64_t)(1019+b)<<52)|((uint64_t)(2*k+1)<<46))
// Bands cover ratio in [1/16, 16], each band spans one octave, 32 cells/band.
static const Atan2LogData atan2_log_table[256] = {
    // b=0: c in [1/16, 1/8)
    {0x1.03a6d1c06693dp-4, -0x1.440a69a6af54cp-131},  // k=0  c=0.015625
    {0x1.0b9e589974c5ap-4, -0x1.eba31cfa47c34p-132},  // k=1  c=0.032227
    {0x1.13955a967a682p-4, -0x1.a21e8f9a94f0ap-130},  // k=2  c=0.033203
    {0x1.1b8bd3d2450b4p-4, -0x1.98144b78c45e3p-133},  // k=3  c=0.034180
    {0x1.2381c06936f53p-4,  0x1.82304a7e47ebcp-133},  // k=4  c=0.035156
    {0x1.2b771c79524f9p-4,  0x1.58a0099e4372fp-130},  // k=5  c=0.036133
    {0x1.336be4224448fp-4,  0x1.ed6510c2c649ep-130},  // k=6  c=0.037109
    {0x1.3b6013857029ap-4,  0x1.15939ecbe466dp-134},  // k=7  c=0.038086
    {0x1.4353a6c5fa5c9p-4, -0x1.1f1dcf4bff49cp-130},  // k=8  c=0.039063
    {0x1.4b469a08d36a9p-4, -0x1.cc226adce4cfcp-131},  // k=9  c=0.040039
    {0x1.5338e974c2e93p-4,  0x1.65414a8fd4b1bp-133},  // k=10 c=0.041016
    {0x1.5b2a9132725bep-4,  0x1.e5376d20e71c8p-130},  // k=11 c=0.041992
    {0x1.631b8d6c78073p-4, -0x1.af41d321d779bp-130},  // k=12 c=0.042969
    {0x1.6b0bda4f61b64p-4, -0x1.4b1e9de48ad55p-130},  // k=13 c=0.043945
    {0x1.72fb7409bf71ep-4,  0x1.5644df963b5b6p-130},  // k=14 c=0.044922
    {0x1.7aea56cc2e292p-4, -0x1.a088e32c2c100p-130},  // k=15 c=0.045898
    {0x1.82d87ec9624b0p-4,  0x1.5adb7cac8503fp-130},  // k=16 c=0.046875
    {0x1.8ac5e8363250cp-4, -0x1.a2590bdf8317fp-132},  // k=17 c=0.047852
    {0x1.92b28f49a1396p-4, -0x1.2ce5845482901p-130},  // k=18 c=0.048828
    {0x1.9a9e703ce8f4bp-4, -0x1.5635b14f37529p-130},  // k=19 c=0.049805
    {0x1.a289874b84bf2p-4,  0x1.ba1fae2e4a314p-131},  // k=20 c=0.050781
    {0x1.aa73d0b33b6cdp-4, -0x1.09d47cb9c124dp-132},  // k=21 c=0.051758
    {0x1.b25d48b429a49p-4, -0x1.e7b8c800b49c9p-130},  // k=22 c=0.052734
    {0x1.ba45eb90cc09dp-4,  0x1.9af102d953d26p-131},  // k=23 c=0.053711
    {0x1.c22db58e0955ep-4, -0x1.a22968a463a9bp-130},  // k=24 c=0.054688
    {0x1.ca14a2f33c5fcp-4, -0x1.7b5d34dccfba7p-131},  // k=25 c=0.055664
    {0x1.d1fab00a3e127p-4,  0x1.2db03ba8429f3p-133},  // k=26 c=0.056641
    {0x1.d9dfd91f6f51fp-4, -0x1.610a1c2c28ef2p-132},  // k=27 c=0.057617
    {0x1.e1c41a81c2cd8p-4, -0x1.f8c949711d1b0p-132},  // k=28 c=0.058594
    {0x1.e9a77082c6c06p-4,  0x1.dc7a2de838284p-130},  // k=29 c=0.059570
    {0x1.f189d776aea02p-4,  0x1.a8acda409f0bcp-130},  // k=30 c=0.060547
    {0x1.f96b4bb45cb78p-4,  0x1.3200206532dfcp-131},  // k=31 c=0.061523
    // b=1: c in [1/8, 1/4)
    {0x1.029dd57ffc493p-3, -0x1.ed8a26a3fa1dcp-132},  // k=0  c=0.125977
    {0x1.0a7c5b4bed20fp-3,  0x1.3409971abb833p-132},  // k=1  c=0.132813
    {0x1.1258daff330b4p-3, -0x1.a02ca65ec4da8p-133},  // k=2  c=0.140625
    {0x1.1a334638df0d3p-3, -0x1.9d9e4e05e3129p-130},  // k=3  c=0.148438
    {0x1.220b8eafa5aa3p-3,  0x1.58135f665d605p-131},  // k=4  c=0.156250
    {0x1.29e1a6326d7d6p-3, -0x1.1a14879f2918cp-131},  // k=5  c=0.164063
    {0x1.31b57ea8db38dp-3,  0x1.02a5c6ee64ffep-130},  // k=6  c=0.171875
    {0x1.39870a13dafd5p-3, -0x1.6b35b3e958d6cp-132},  // k=7  c=0.179688
    {0x1.41563a8e2700dp-3,  0x1.6d833c9d54ea7p-130},  // k=8  c=0.187500
    {0x1.4923024ccb781p-3,  0x1.e9e3ee112bae1p-130},  // k=9  c=0.195313
    {0x1.50ed539fa7b92p-3, -0x1.e9e3cf573184cp-132},  // k=10 c=0.203125
    {0x1.58b520f1ec8e1p-3,  0x1.140c1041192b2p-133},  // k=11 c=0.210938
    {0x1.607a5cca97ad8p-3, -0x1.e073fabff9aa8p-131},  // k=12 c=0.218750
    {0x1.683cf9ccec514p-3,  0x1.865450cc24ef2p-132},  // k=13 c=0.226563
    {0x1.6ffceab8e8e2cp-3,  0x1.1744c2103a78dp-131},  // k=14 c=0.234375
    {0x1.77ba226bb9b5ap-3, -0x1.29531bc9ac913p-131},  // k=15 c=0.242188
    {0x1.7f7493e028c98p-3,  0x1.ca1ace3a65c63p-130},  // k=16 c=0.250000
    {0x1.872c322f0a8ccp-3, -0x1.1cfe2a9ea0444p-130},  // k=17 c=0.257813
    {0x1.8ee0f08fa79a2p-3,  0x1.8a2a50dd0175ap-131},  // k=18 c=0.265625
    {0x1.9692c258236b8p-3, -0x1.804bd533c09d6p-131},  // k=19 c=0.273438
    {0x1.9e419afddffe1p-3,  0x1.7f7930120d408p-130},  // k=20 c=0.281250
    {0x1.a5ed6e15de61fp-3,  0x1.5e889b1134f36p-131},  // k=21 c=0.289063
    {0x1.ad962f551c32fp-3,  0x1.79780d3bc7200p-130},  // k=22 c=0.296875
    {0x1.b53bd290edf69p-3,  0x1.295e234391dbcp-130},  // k=23 c=0.304688
    {0x1.bcde4bbf565c2p-3,  0x1.daa0128157d4bp-130},  // k=24 c=0.312500
    {0x1.c47d8ef75a5dcp-3,  0x1.993b1d3f03482p-130},  // k=25 c=0.320313
    {0x1.cc199071523f5p-3, -0x1.28dfdf79206b8p-131},  // k=26 c=0.328125
    {0x1.d3b24487376b7p-3, -0x1.e9abfda455927p-130},  // k=27 c=0.335938
    {0x1.db479fb4ef2c9p-3,  0x1.2c48444847e72p-130},  // k=28 c=0.343750
    {0x1.e2d996989242ep-3, -0x1.9501ca5ac4181p-130},  // k=29 c=0.351563
    {0x1.ea681df2b156bp-3,  0x1.0117938736292p-130},  // k=30 c=0.359375
    {0x1.f1f32aa696486p-3,  0x1.6c2d7211f25d5p-130},  // k=31 c=0.367188
    // b=2: c in [1/4, 1/2)
    {0x1.fd3d1fc40dbe4p-3,  0x1.a27bf50da7f87p-131},  // k=0  c=0.257813
    {0x1.061eea03d6291p-2, -0x1.c044a1c1bd879p-137},  // k=1  c=0.265625
    {0x1.0d97ee509acb3p-2, -0x1.94300afd7a8bbp-130},  // k=2  c=0.273438
    {0x1.150973a9ce547p-2,  0x1.2fc5c9b3fd080p-130},  // k=3  c=0.281250
    {0x1.1c735212dd884p-2,  0x1.131740ef52becp-132},  // k=4  c=0.289063
    {0x1.23d562b381042p-2, -0x1.7cbdbcc2e3e29p-130},  // k=5  c=0.296875
    {0x1.2b2f7fd9b5fe2p-2, -0x1.661cb5aa5d4e2p-133},  // k=6  c=0.304688
    {0x1.328184fb58952p-2, -0x1.7820fae8a08afp-130},  // k=7  c=0.312500
    {0x1.39cb4eb76157cp-2,  0x1.aa6b1c0e96258p-132},  // k=8  c=0.320313
    {0x1.410cbad6c7d33p-2, -0x1.94b47d40dbf0fp-131},  // k=9  c=0.328125
    {0x1.4845a84d0c21bp-2, -0x1.a75948ace2e7ap-131},  // k=10 c=0.335938
    {0x1.4f75f73869979p-2,  0x1.61f133d0c0908p-130},  // k=11 c=0.343750
    {0x1.569d88e1b4cd8p-2,  0x1.76ee20b4cf8f4p-130},  // k=12 c=0.351563
    {0x1.5dbc3fbbe768dp-2,  0x1.65ce248742f67p-131},  // k=13 c=0.359375
    {0x1.64d1ff635c1c6p-2, -0x1.edd41a72076c1p-130},  // k=14 c=0.367188
    {0x1.6bdeac9cbd76dp-2,  0x1.5da512b28ee3dp-130},  // k=15 c=0.375000
    {0x1.72e22d53aa2aap-2,  0x1.efa88b6e05513p-131},  // k=16 c=0.382813
    {0x1.79dc6899118d1p-2, -0x1.d527c3b8ae77cp-132},  // k=17 c=0.390625
    {0x1.80cd46a14b1d1p-2,  0x1.8630a3e0d691fp-130},  // k=18 c=0.398438
    {0x1.87b4b0c1ebedcp-2, -0x1.53fd4584528c1p-131},  // k=19 c=0.406250
    {0x1.8e92916f5cde8p-2, -0x1.1f22996159aa8p-132},  // k=20 c=0.414063
    {0x1.9566d43a34907p-2, -0x1.a307ad86628bcp-130},  // k=21 c=0.421875
    {0x1.9c3165cc58107p-2,  0x1.4e08393f865a3p-131},  // k=22 c=0.429688
    {0x1.a2f233e5e530bp-2,  0x1.9ca2064123d21p-130},  // k=23 c=0.437500
    {0x1.a9a92d59e98cfp-2, -0x1.e0371529fefb9p-131},  // k=24 c=0.445313
    {0x1.b056420ae9344p-2, -0x1.94aa1b63c8096p-137},  // k=25 c=0.453125
    {0x1.b6f962e737efcp-2, -0x1.9435b03cc4d35p-131},  // k=26 c=0.460938
    {0x1.bd9281e528192p-2, -0x1.92d084ba29ab2p-130},  // k=27 c=0.468750
    {0x1.c42191ff11eb7p-2,  0x1.a36e707bdd985p-130},  // k=28 c=0.476563
    {0x1.caa6872f3631bp-2,  0x1.5864e6166008ep-130},  // k=29 c=0.484375
    {0x1.d121566b7f2adp-2, -0x1.8a2604190442bp-131},  // k=30 c=0.492188
    {0x1.d791f5a1226f5p-2, -0x1.a5d259690d838p-131},  // k=31 c=0.500000
    // b=3: c in [1/2, 1)
    {0x1.e127b6b0744b0p-2, -0x1.c3c5d4084ef34p-131},  // k=0  c=0.515625
    {0x1.edcb6d43f8435p-2,  0x1.6042ea30dc440p-131},  // k=1  c=0.531250
    {0x1.fa45dd3029259p-2, -0x1.94b09934ab218p-132},  // k=2  c=0.546875
    {0x1.034b709250488p-1, -0x1.b5c586fcdef8ep-131},  // k=3  c=0.562500
    {0x1.095f30861a590p-1, -0x1.456e37be8d0a3p-130},  // k=4  c=0.578125
    {0x1.0f5e28b67e295p-1, -0x1.908a50183a518p-130},  // k=5  c=0.593750
    {0x1.154859637646ap-1, -0x1.89081332152abp-136},  // k=6  c=0.609375
    {0x1.1b1dc87904285p-1, -0x1.01ea3cf51279fp-131},  // k=7  c=0.625000
    {0x1.20de813e823b2p-1, -0x1.983d6deee1590p-130},  // k=8  c=0.640625
    {0x1.268a940696da6p-1, -0x1.d9258369d58a1p-131},  // k=9  c=0.656250
    {0x1.2c2215e024466p-1,  0x1.26bcee46281a9p-130},  // k=10 c=0.671875
    {0x1.31a52048874bep-1, -0x1.19f7b82d31fcfp-135},  // k=11 c=0.687500
    {0x1.3713d0df6c504p-1, -0x1.d12a9f8932550p-131},  // k=12 c=0.703125
    {0x1.3c6e491c78dc5p-1,  0x1.ca7bd64272a9cp-131},  // k=13 c=0.718750
    {0x1.41b4ae06fea41p-1,  0x1.cb31765f46ca6p-131},  // k=14 c=0.734375
    {0x1.46e727efe4716p-1, -0x1.f221784eca08bp-130},  // k=15 c=0.750000
    {0x1.4c05e22de94e5p-1,  0x1.190993172408fp-130},  // k=16 c=0.765625
    {0x1.51110adc5ed81p-1, -0x1.e854b5a46ffd6p-132},  // k=17 c=0.781250
    {0x1.5608d29c70c34p-1, -0x1.7891861e53bf3p-134},  // k=18 c=0.796875
    {0x1.5aed6c5909517p-1,  0x1.4a7d3d0d55fbap-130},  // k=19 c=0.812500
    {0x1.5fbf0d0d5cc4ap-1,  0x1.f9044c679a55bp-130},  // k=20 c=0.828125
    {0x1.647deb8e20b90p-1, -0x1.357108d67a7d4p-130},  // k=21 c=0.843750
    {0x1.692a40556fb6ap-1,  0x1.d09bceee4046ap-130},  // k=22 c=0.859375
    {0x1.6dc44551553afp-1,  0x1.5af43d90c57d5p-133},  // k=23 c=0.875000
    {0x1.724c35b4fae7bp-1,  0x1.7f7a0cffe29bap-134},  // k=24 c=0.890625
    {0x1.76c24dcc6c6c0p-1,  0x1.342deff60661ep-132},  // k=25 c=0.906250
    {0x1.7b26cad2e50fep-1, -0x1.d8138764fd9cdp-130},  // k=26 c=0.921875
    {0x1.7f79eacb97898p-1, -0x1.8efbd105c6a7fp-131},  // k=27 c=0.937500
    {0x1.83bbec5cdee22p-1,  0x1.812b904f1f352p-132},  // k=28 c=0.953125
    {0x1.87ed0eadc5a2ap-1, -0x1.8daacee1c531ep-131},  // k=29 c=0.968750
    {0x1.8c0d9145cf49dp-1,  0x1.bdd2d76a28129p-130},  // k=30 c=0.984375
    {0x1.901db3eeef187p-1,  0x1.dbe94df613465p-131},  // k=31 c=1.000000
    // b=4: c in [1, 2)
    {0x1.9617bfeebaf28p-1,  0x1.0032c8866f9b2p-56},   // k=0  c=1.01563
    {0x1.9dd8d51585769p-1, -0x1.ae877ffec0468p-55},   // k=1  c=1.04688
    {0x1.a55ce84fecd16p-1, -0x1.21ee7dfa97e23p-55},   // k=2  c=1.07813
    {0x1.aca5f3f19ab96p-1,  0x1.e8e23fc4ff624p-55},   // k=3  c=1.10938
    {0x1.b3b5ec3794397p-1, -0x1.2985d9907b080p-55},   // k=4  c=1.14063
    {0x1.ba8ebd7c5a191p-1,  0x1.e58aede56b11ep-55},   // k=5  c=1.17188
    {0x1.c1324abe8915dp-1,  0x1.5de7930719b20p-61},   // k=6  c=1.20313
    {0x1.c7a26c706c612p-1, -0x1.09a8c1f8a4817p-56},   // k=7  c=1.23438
    {0x1.cde0ef87a615dp-1,  0x1.88b5df85eda52p-56},   // k=8  c=1.26563
    {0x1.d3ef94c4cabf9p-1,  0x1.ea50304a31a71p-57},   // k=9  c=1.29688
    {0x1.d9d0102b87e53p-1,  0x1.9dff587921d15p-55},   // k=10 c=1.32813
    {0x1.df8408a4a4bf5p-1,  0x1.0f428b385436fp-55},   // k=11 c=1.35938
    {0x1.e50d17c3dc497p-1, -0x1.844a35b61a4b4p-55},   // k=12 c=1.39063
    {0x1.ea6cc9ac363efp-1,  0x1.64f83238896e5p-58},   // k=13 c=1.42188
    {0x1.efa49d0e1fcefp-1,  0x1.85c3eaf843caap-55},   // k=14 c=1.45313
    {0x1.f4b6033b14105p-1, -0x1.c5354f363c7eep-58},   // k=15 c=1.48438
    {0x1.f9a2604b273d7p-1,  0x1.59fb0e0020faap-55},   // k=16 c=1.51563
    {0x1.fe6b0b513e96ap-1, -0x1.82c681ffe40a2p-56},   // k=17 c=1.54688
    {0x1.0188a74d94efcp+0, -0x1.de8df32afb52dp-54},   // k=18 c=1.57813
    {0x1.03cb33fd99ad1p+0,  0x1.97daa58fd9fa9p-54},   // k=19 c=1.60938
    {0x1.05fdc48c884c5p+0,  0x1.1c07d0e3c5400p-54},   // k=20 c=1.64063
    {0x1.0820ebe435222p+0, -0x1.f36f8ddb2aea1p-54},   // k=21 c=1.67188
    {0x1.0a353731b9f23p+0,  0x1.95b45bbe0eb51p-55},   // k=22 c=1.70313
    {0x1.0c3b2e1423c0fp+0, -0x1.ac7ddfb7312b5p-55},   // k=23 c=1.73438
    {0x1.0e3352cbf6450p+0, -0x1.4e55e38fb33b0p-54},   // k=24 c=1.76563
    {0x1.101e226b0f2b4p+0, -0x1.ca14e83b1fd9ap-55},   // k=25 c=1.79688
    {0x1.11fc150486eeap+0,  0x1.b640c9c014a5ap-54},   // k=26 c=1.82813
    {0x1.13cd9ddc3dfe6p+0,  0x1.00949f8187484p-56},   // k=27 c=1.85938
    {0x1.15932b95d395fp+0, -0x1.ba87ea99b85f1p-54},   // k=28 c=1.89063
    {0x1.174d2862d04ccp+0,  0x1.dff7f9a86ea82p-54},   // k=29 c=1.92188
    {0x1.18fbfa2fd93b8p+0,  0x1.6fe2c91a1ae8dp-55},   // k=30 c=1.95313
    {0x1.1aa002d0c9dabp+0,  0x1.b3a673c04199ap-55},   // k=31 c=1.98438
    // b=5: c in [2, 4)
    {0x1.1d02a2edce14fp+0,  0x1.d0a3b2bae2099p-56},   // k=0  c=2.03125
    {0x1.200e5ae0dd61dp+0,  0x1.bed3cdd901e18p-55},   // k=1  c=2.09375
    {0x1.22f5127ff70b5p+0,  0x1.4a430e2d0a41ep-55},   // k=2  c=2.15625
    {0x1.25b92ee2f7a49p+0, -0x1.75e638d53460ap-54},   // k=3  c=2.21875
    {0x1.285ce59c97320p+0, -0x1.75b26f6f4d78fp-54},   // k=4  c=2.28125
    {0x1.2ae240bd32395p+0, -0x1.fd6bc679da199p-54},   // k=5  c=2.34375
    {0x1.2d4b22823a35cp+0,  0x1.cf7b124d9fd08p-56},   // k=6  c=2.40625
    {0x1.2f9948b6a1b4cp+0,  0x1.22febbf47b080p-55},   // k=7  c=2.46875
    {0x1.31ce4fc931347p+0,  0x1.39a6d04c1b60bp-54},   // k=8  c=2.53125
    {0x1.33ebb59df782bp+0, -0x1.02b3d3b7bf8dcp-54},   // k=9  c=2.59375
    {0x1.35f2dc20141bbp+0, -0x1.f41f13d3ea682p-54},   // k=10 c=2.65625
    {0x1.37e50b98fe91cp+0,  0x1.7c4f65ff85e9ep-54},   // k=11 c=2.71875
    {0x1.39c374d238f7bp+0, -0x1.0e10bce7c08c4p-57},   // k=12 c=2.78125
    {0x1.3b8f3306167a9p+0, -0x1.15725c9547d04p-54},   // k=13 c=2.84375
    {0x1.3d494da3f1dcap+0, -0x1.39e5cc1bf2c1bp-54},   // k=14 c=2.90625
    {0x1.3ef2b9ebde493p+0,  0x1.0ab7d20ab732ap-54},   // k=15 c=2.96875
    {0x1.408c5c658b763p+0, -0x1.c066d31f57e75p-55},   // k=16 c=3.03125
    {0x1.42170a35c681ap+0,  0x1.e5c774907e6dbp-55},   // k=17 c=3.09375
    {0x1.43938a55b4f1dp+0,  0x1.1ff9682e290acp-54},   // k=18 c=3.15625
    {0x1.450296ae9f949p+0, -0x1.cfd6fc36c6974p-54},   // k=19 c=3.21875
    {0x1.4664dd1ce02f9p+0,  0x1.38ecc63a07d2dp-54},   // k=20 c=3.28125
    {0x1.47bb005c4732dp+0,  0x1.a2f98d0b4820dp-54},   // k=21 c=3.34375
    {0x1.490598e015e2ep+0,  0x1.e07def52831b8p-54},   // k=22 c=3.40625
    {0x1.4a453598759c0p+0, -0x1.666bc3f7a196ep-55},   // k=23 c=3.46875
    {0x1.4b7a5ca7259aap+0, -0x1.ec80f282dc460p-54},   // k=24 c=3.53125
    {0x1.4ca58c04ede32p+0,  0x1.144de97a969fcp-55},   // k=25 c=3.59375
    {0x1.4dc73a193f1b4p+0,  0x1.12518d76d3566p-54},   // k=26 c=3.65625
    {0x1.4edfd645441e7p+0, -0x1.f67eaeec21278p-54},   // k=27 c=3.71875
    {0x1.4fefc9638a79fp+0,  0x1.ca0501d485c89p-54},   // k=28 c=3.78125
    {0x1.50f7763d4a801p+0, -0x1.d0bda250c54e3p-54},   // k=29 c=3.84375
    {0x1.51f739f63df6ap+0, -0x1.34fffe9d872fbp-54},   // k=30 c=3.90625
    {0x1.52ef6c6fdd517p+0, -0x1.ef7a6346bb398p-54},   // k=31 c=3.96875
    // b=6: c in [4, 8)
    {0x1.54563c118794cp+0,  0x1.dcf8852691640p-54},   // k=0  c=4.06250
    {0x1.561d0ad75b853p+0, -0x1.3c563b0cd09ecp-61},   // k=1  c=4.18750
    {0x1.57cade57dba51p+0, -0x1.8b8666c0ed546p-60},   // k=2  c=4.31250
    {0x1.5961aeac184b3p+0, -0x1.ae4b95d52e176p-54},   // k=3  c=4.43750
    {0x1.5ae3412f1c467p+0,  0x1.41f255fc8a11ap-54},   // k=4  c=4.56250
    {0x1.5c512e9c5bb9fp+0, -0x1.cf8fbd02712bep-54},   // k=5  c=4.68750
    {0x1.5dace854dff1cp+0,  0x1.73f5a5a50fadfp-59},   // k=6  c=4.81250
    {0x1.5ef7bcecde57ap+0,  0x1.c1fc96bde482ap-54},   // k=7  c=4.93750
    {0x1.6032dc1db31f6p+0, -0x1.90a89eba7946cp-55},   // k=8  c=5.06250
    {0x1.615f5a338cc04p+0,  0x1.a791fa40522aap-54},   // k=9  c=5.18750
    {0x1.627e330a40ae9p+0,  0x1.2a8b12af9f40bp-58},   // k=10 c=5.31250
    {0x1.63904ca99b192p+0,  0x1.0637e25b5a21fp-57},   // k=11 c=5.43750
    {0x1.6496798ee1ec5p+0, -0x1.62b32e1c7b9dbp-56},   // k=12 c=5.56250
    {0x1.65917aaf1caa2p+0, -0x1.90c4cafa9df20p-54},   // k=13 c=5.68750
    {0x1.6682013be9a64p+0, -0x1.9479fbb399bacp-58},   // k=14 c=5.81250
    {0x1.6768b0332cd22p+0, -0x1.053d6aa6ae0b1p-54},   // k=15 c=5.93750
    {0x1.68461dc0a6bc7p+0, -0x1.d5023b5c99699p-54},   // k=16 c=6.06250
    {0x1.691ad47778958p+0,  0x1.6446b89f76d5bp-54},   // k=17 c=6.18750
    {0x1.69e75468badcdp+0,  0x1.f296eed20e524p-55},   // k=18 c=6.31250
    {0x1.6aac141b90791p+0, -0x1.ae605e48c3e5ap-54},   // k=19 c=6.43750
    {0x1.6b69816a8160fp+0, -0x1.acc41f3aad1cbp-55},   // k=20 c=6.56250
    {0x1.6c20024961b81p+0,  0x1.f035f98c5e92dp-54},   // k=21 c=6.68750
    {0x1.6ccff57698526p+0,  0x1.495b739d3b0cfp-56},   // k=22 c=6.81250
    {0x1.6d79b31a369dfp+0, -0x1.81895cac911d7p-54},   // k=23 c=6.93750
    {0x1.6e1d8d5501419p+0,  0x1.cb5ebec61d812p-54},   // k=24 c=7.06250
    {0x1.6ebbd0c142077p+0,  0x1.1366b90f3db0bp-55},   // k=25 c=7.18750
    {0x1.6f54c4e6ff0e5p+0, -0x1.2c3eb2e5d4742p-54},   // k=26 c=7.31250
    {0x1.6fe8aca4ff2b9p+0,  0x1.aca07cf67f4ccp-54},   // k=27 c=7.43750
    {0x1.7077c68fd68fcp+0, -0x1.0b6cc09dfda37p-55},   // k=28 c=7.56250
    {0x1.71024d48100cfp+0, -0x1.29dfd27c252d8p-54},   // k=29 c=7.68750
    {0x1.718877c865e8ap+0,  0x1.668bdc523dcecp-54},   // k=30 c=7.81250
    {0x1.720a79ace01eep+0,  0x1.969999e452063p-54},   // k=31 c=7.93750
    // b=7: c in [8, 16)
    {0x1.72c619a79eb51p+0, -0x1.f7e5d0ba64e0dp-54},   // k=0  c=8.12500
    {0x1.73b36322e5a78p+0, -0x1.b27f69dd145cap-54},   // k=1  c=8.37500
    {0x1.74931b40e6089p+0,  0x1.3b1ccbff02344p-54},   // k=2  c=8.62500
    {0x1.7566626754235p+0, -0x1.79b006f7943c1p-54},   // k=3  c=8.87500
    {0x1.762e3a0ad27dbp+0,  0x1.0fce0b5e61c3ep-54},   // k=4  c=9.12500
    {0x1.76eb88b56e0e1p+0, -0x1.6e6b9d43c04eap-56},   // k=5  c=9.37500
    {0x1.779f1d70e8bfdp+0, -0x1.282f0950761e8p-55},   // k=6  c=9.62500
    {0x1.7849b2afa9392p+0,  0x1.03bfbbf57e392p-54},   // k=7  c=9.87500
    {0x1.78ebf0ca06a41p+0,  0x1.a55484f48da49p-54},   // k=8  c=10.1250
    {0x1.798670219969dp+0,  0x1.78e1236d13372p-56},   // k=9  c=10.3750
    {0x1.7a19baf900435p+0, -0x1.f45ce761c29b9p-56},   // k=10 c=10.6250
    {0x1.7aa64f0bf4bbep+0,  0x1.7d1bd08e87443p-54},   // k=11 c=10.8750
    {0x1.7b2c9ef177685p+0,  0x1.ce4853ff4e5a8p-58},   // k=12 c=11.1250
    {0x1.7bad13502f261p+0, -0x1.033d85ad876a7p-54},   // k=13 c=11.3750
    {0x1.7c280bebba7a6p+0, -0x1.ef3f029f07962p-54},   // k=14 c=11.6250
    {0x1.7c9de0909622ap+0, -0x1.df66031923dd3p-55},   // k=15 c=11.8750
    {0x1.7d0ee1e3533e4p+0, -0x1.f0838cccc1413p-55},   // k=16 c=12.1250
    {0x1.7d7b5a1718148p+0,  0x1.00655b0b0e0b2p-57},   // k=17 c=12.3750
    {0x1.7de38d8ec8635p+0,  0x1.4d61badd75b84p-55},   // k=18 c=12.6250
    {0x1.7e47bb6baf0f3p+0, -0x1.17f759533be11p-54},   // k=19 c=12.8750
    {0x1.7ea81e0c15f4ap+0,  0x1.7d7da6fa30a02p-54},   // k=20 c=13.1250
    {0x1.7f04eb7bdd7d0p+0,  0x1.4ec9fb8a41b13p-54},   // k=21 c=13.3750
    {0x1.7f5e55d8d9377p+0, -0x1.6fd26c21f97a7p-54},   // k=22 c=13.6250
    {0x1.7fb48bac767e9p+0,  0x1.23de9f42ff0c2p-54},   // k=23 c=13.8750
    {0x1.8007b83bfc040p+0, -0x1.4e96acfd0c2aep-57},   // k=24 c=14.1250
    {0x1.805803d0833eap+0, -0x1.527c3fb90a325p-61},   // k=25 c=14.3750
    {0x1.80a593f7a5d00p+0,  0x1.c60b1e17effadp-56},   // k=26 c=14.6250
    {0x1.80f08bbdb7d12p+0, -0x1.4827d93ae77aap-57},   // k=27 c=14.8750
    {0x1.81390be24bb09p+0, -0x1.0577a17c47468p-54},   // k=28 c=15.1250
    {0x1.817f3307a416ep+0, -0x1.71f4e2a83101ap-54},   // k=29 c=15.3750
    {0x1.81c31ddda3809p+0, -0x1.6cc7075443c24p-56},   // k=30 c=15.6250
    {0x1.8204e748b75e2p+0, -0x1.2aa5b7c905744p-54},   // k=31 c=15.8750
};

// Main-path polynomial: atan(r) = r + r*s*(A0 + s*(A1 + s*A2)), s = r*r
// Valid for |r| <= 1/128. Minimax over that interval.
static const double ATAN2_A0 = -0x1.5555555555547p-2;
static const double ATAN2_A1 =  0x1.9999999199c0fp-3;
static const double ATAN2_A2 = -0x1.248cf3e342e15p-3;

// Fallback polynomial: atan(z) = z + z*s*(B0 + s*(B1 + ... + s*B5)), s = z*z
// Valid for |z| <= 1/16. Minimax over that interval.
static const double ATAN2_B0 = -0x1.5555555555555p-2;
static const double ATAN2_B1 =  0x1.999999999994ep-3;
static const double ATAN2_B2 = -0x1.249249245adbfp-3;
static const double ATAN2_B3 =  0x1.c71c6fedda3b7p-4;
static const double ATAN2_B4 = -0x1.74598599ff148p-4;
static const double ATAN2_B5 =  0x1.37e5cbb2db80ep-4;

// Shared table indexing: cells per octave band
#define ATAN2_LOG2SIZE 5
// atan2/atan2f fallback cutoffs: unbiased exponent of ratio
#define ATAN2_FALLBACK_EMIN_DP64 -53  // (z+1)*z == z in double: requires |z| < ulp(1.0)/2 = 2^-53
#define ATAN2_FALLBACK_EMAX_DP64  53  // pi/2 - z == pi/2 in double: requires z < ulp(pi/2)/2 = 2^-53
#define ATAN2_FALLBACK_EMIN_SP32 -53  // (z+1)*z == z in double: requires |z| < ulp(1.0)/2 = 2^-53
#define ATAN2_FALLBACK_EMAX_SP32  25  // pi/2 - z == pi/2 in float: requires z < ulp(pi/2_f)/2 = 2^-24, i.e. z <= 2^-25

// Double-precision pi constants shared by atan2 and atan2f
static const double atan2_pi      = 0x1.921fb54442d18p+1;
static const double atan2_piby2   = 0x1.921fb54442d18p+0;
static const double atan2_piby4   = 0x1.921fb54442d18p-1;
static const double atan2_pi3by4  = 0x1.2d97c7f3321d2p+1;
static const double atan2_pi_head = 0x1.921fb50000000p+1;
static const double atan2_pi_tail = 0x1.110b4611a6263p-25;  // pi - pi_head

#endif /* ATAN2_DATA_H */
