//
// Created by Rain Night on 2024/1/11.
//
#include "org_apache_hadoop_hbase_util_RibbonHelper.h"
using TypeParam = Settings_Coeff128_Homog;
IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
IMPORT_RIBBON_IMPL_TYPES(TypeParam);

// 全局变量声明
InterleavedSoln *isoln = nullptr;
Banding *banding = nullptr;
Hasher *hasher = nullptr;

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    initRibbonFilter
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_initRibbonFilter
        (JNIEnv *, jobject, jint ji) {
//    Index num_slots = 131072;
    Index num_slots = static_cast<uint32_t>(ji);
//    uint32_t num_to_add = 109306;
    uint32_t max_ibytes = static_cast<uint32_t>(sizeof(ResultRow) * num_slots);
    std::unique_ptr<char[]> idata(new char[max_ibytes]);
    isoln = new InterleavedSoln(idata.get(), max_ibytes);
    banding = new Banding();
    banding->SetOrdinalSeed(0);
    banding->Reset(num_slots);
    hasher = new Hasher();
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    addKey
 * Signature: (Ljava/lang/String;)Z
 */
JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_addKey
        (JNIEnv * env, jobject jo, jstring js) {
    const char *cString = env->GetStringUTFChars(js , nullptr);
    return banding->Add(Slice(cString));
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    backSubst
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_backSubst
        (JNIEnv *, jobject) {
    isoln->BackSubstFrom(*banding);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    filterQuery
 * Signature: (Ljava/lang/String;)Z
 */
JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_filterQuery
        (JNIEnv * env, jobject jo, jstring js) {
    const char *cString = env->GetStringUTFChars(js , nullptr);
    return isoln->FilterQuery(Slice(cString), *hasher);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    close
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_close
        (JNIEnv *, jobject) {
    delete banding;
    delete isoln;
    delete hasher;
}